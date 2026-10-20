"""Deploy Prefect flows from a config file.

Every deployment is described in gb_prefect/deployments/deployments.yaml. Your own paths 
in deployments.local.yaml next to it, or are passed with --set name=value.

The worker clones --branch from GitHub at run time, so commit and push the branch first.

The whole config is resolved and validated before anything is deployed: a missing variable,
an unknown key, pool or deployment name stops the run with nothing deployed.

Usage: run from the repository root, with PREFECT_API_URL pointing at the Prefect server

    # show what would be deployed, deploy nothing
    python gb_prefect/deployments/deploy.py --branch <branch> --dry-run

    # deploy everything
    python gb_prefect/deployments/deploy.py --branch <branch>

    # deploy two of them, as test copies, with a different output directory
    python gb_prefect/deployments/deploy.py --branch <branch> \\
        --only busco-gca busco-dispatch-automatic --name-suffix _test \\
        --set busco_outdir=/path/to/busco_output
"""

import argparse
import copy
import json
import re
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml  # type: ignore

DEPLOYMENTS_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG = DEPLOYMENTS_DIR / "deployments.yaml"
DEFAULT_LOCAL_CONFIG = DEPLOYMENTS_DIR / "deployments.local.yaml"

TOP_LEVEL_KEYS = {"repo", "variables", "work_pools", "deployments"}
POOL_KEYS = {"job_variables"}
DEPLOYMENT_KEYS = {
    "entrypoint",
    "work_pool",
    "parameters",
    "job_variables",
    "tags",
    "description",
    "concurrency_limit",
}
REQUIRED_DEPLOYMENT_KEYS = {"entrypoint", "work_pool"}
BASE_TAG = "genebuild"

VARIABLE_PATTERN = re.compile(r"\{(\w+)\}")


class ConfigError(ValueError):
    """The deployment config is invalid; the message lists every problem found."""


def deep_merge(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    """Return base with override merged in: dicts merged key by key, other values replaced."""
    merged = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = deep_merge(merged[key], value)
        else:
            merged[key] = copy.deepcopy(value)
    return merged


def load_config(config_path: Path, local_path: Optional[Path]) -> Dict[str, Any]:
    """Load config_path, with local_path (if given and present) merged on top."""
    config = yaml.safe_load(config_path.read_text()) or {}
    if local_path is not None and local_path.exists():
        config = deep_merge(config, yaml.safe_load(local_path.read_text()) or {})
    return config


def parse_set_args(pairs: List[str]) -> Dict[str, str]:
    """Turn ["name=value", ...] from --set into a dict."""
    variables = {}
    for pair in pairs:
        name, sep, value = pair.partition("=")
        if not sep or not name:
            raise ConfigError(f"--set expects name=value, got '{pair}'")
        variables[name] = value
    return variables


def substitute(value: Any, variables: Dict[str, Any], missing: set) -> Any:
    """Replace "{name}" in every string inside value; names not set are added to missing."""
    if isinstance(value, str):

        def replace(match: "re.Match[str]") -> str:
            name = match.group(1)
            if variables.get(name) is None:
                missing.add(name)
                return match.group(0)
            return str(variables[name])

        return VARIABLE_PATTERN.sub(replace, value)
    if isinstance(value, dict):
        return {key: substitute(item, variables, missing) for key, item in value.items()}
    if isinstance(value, list):
        return [substitute(item, variables, missing) for item in value]
    return value


def _unknown_keys(section: str, data: Any, allowed: set) -> List[str]:
    if not isinstance(data, dict):
        return [f"{section} must be a mapping"]
    return [f"{section}: unknown key '{key}'" for key in sorted(set(data) - allowed)]


def build_specs(  # pylint: disable=too-many-locals
    config: Dict[str, Any],
    branch: str,
    name_suffix: str = "",
    only: Optional[List[str]] = None,
) -> List[Dict[str, Any]]:
    """Resolve config into one deploy() spec per selected deployment.

    Raises ConfigError listing every problem (not just the first) if anything is wrong, so
    a deploy never stops halfway through.
    """
    errors = _unknown_keys("config", config, TOP_LEVEL_KEYS)
    variables = config.get("variables") or {}
    pools = config.get("work_pools") or {}
    deployments = config.get("deployments") or {}
    if not config.get("repo"):
        errors.append("config: 'repo' is required")
    for pool_name, pool in pools.items():
        errors += _unknown_keys(f"work_pools.{pool_name}", pool or {}, POOL_KEYS)

    selected = only or list(deployments)
    errors += [f"--only: unknown deployment '{name}'" for name in selected if name not in deployments]

    specs = []
    missing: set = set()
    for name in selected:
        deployment = deployments.get(name)
        if deployment is None:
            continue
        section = f"deployments.{name}"
        errors += _unknown_keys(section, deployment, DEPLOYMENT_KEYS)
        errors += [
            f"{section}: '{key}' is required" for key in sorted(REQUIRED_DEPLOYMENT_KEYS - set(deployment))
        ]
        pool_name = deployment.get("work_pool")
        if pool_name not in pools:
            errors.append(f"{section}: work pool '{pool_name}' is not listed under work_pools")
            continue

        job_variables = {
            **((pools[pool_name] or {}).get("job_variables") or {}),
            **(deployment.get("job_variables") or {}),
        }
        description = deployment.get("description", name)
        spec = {
            "name": f"{name}{name_suffix}",
            "entrypoint": deployment.get("entrypoint"),
            "work_pool_name": pool_name,
            "parameters": deployment.get("parameters") or None,
            "job_variables": job_variables or None,
            "tags": [BASE_TAG, *(deployment.get("tags") or []), branch],
            "description": f"{description} from branch {branch}.",
        }
        if deployment.get("concurrency_limit") is not None:
            spec["concurrency_limit"] = deployment["concurrency_limit"]
        specs.append(substitute(spec, variables, missing))

    errors += [
        f"variable '{name}' is not set -- set it in deployments.local.yaml or with --set {name}=..."
        for name in sorted(missing)
    ]
    if errors:
        raise ConfigError("Invalid deployment config:\n  " + "\n  ".join(errors))
    return specs


def warn_if_tracked(path: Path) -> None:
    """Warn if the local override file (meant to stay out of git) is tracked by git."""
    result = subprocess.run(
        ["git", "ls-files", "--error-unmatch", path.name],
        cwd=path.parent,
        capture_output=True,
        check=False,
    )
    if result.returncode == 0:
        print(f"WARNING: {path} is tracked by git -- it is meant to hold local paths only.")


def deploy(specs: List[Dict[str, Any]], repo: str, branch: str) -> None:
    """Create or update each deployment on the Prefect server (PREFECT_API_URL)."""
    # pylint: disable=import-outside-toplevel
    from prefect import flow  # type: ignore
    from prefect.runner.storage import GitRepository  # type: ignore

    source = GitRepository(url=repo, branch=branch)
    for spec in specs:
        options = {key: value for key, value in spec.items() if key != "entrypoint"}
        deployment_id = flow.from_source(source=source, entrypoint=spec["entrypoint"]).deploy(**options)
        print(f"Deployed {spec['entrypoint']} as {spec['name']} ({deployment_id})")


def main() -> None:
    """Parse arguments, resolve the config and deploy (or print it with --dry-run)."""
    parser = argparse.ArgumentParser(description="Deploy Prefect flows from a config file.")
    parser.add_argument("--branch", required=True, help="Git branch the worker clones.")
    parser.add_argument(
        "--only", nargs="+", metavar="NAME", help="Deploy only these deployments (default: all)."
    )
    parser.add_argument(
        "--name-suffix",
        default="",
        help="Appended to every deployment name; leave empty to replace the production ones.",
    )
    parser.add_argument(
        "--set",
        dest="set_vars",
        nargs="+",
        default=[],
        metavar="NAME=VALUE",
        help="Override config variables, e.g. --set outdir=/path enscode=/path.",
    )
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG, help="Config file.")
    parser.add_argument(
        "--local",
        type=Path,
        default=DEFAULT_LOCAL_CONFIG,
        help="Local override file merged on top of --config, if it exists.",
    )
    parser.add_argument(
        "--dry-run", action="store_true", help="Print the resolved deployments; deploy nothing."
    )
    args = parser.parse_args()

    config = load_config(args.config, args.local)
    if args.local.exists():
        warn_if_tracked(args.local)
    try:
        config = deep_merge(config, {"variables": parse_set_args(args.set_vars)})
        specs = build_specs(config, args.branch, args.name_suffix, args.only)
    except ConfigError as err:
        raise SystemExit(str(err)) from err

    if args.dry_run:
        print(json.dumps(specs, indent=2))
        print(f"Dry run: {len(specs)} deployment(s) resolved, nothing deployed.")
        return
    deploy(specs, config["repo"], args.branch)


if __name__ == "__main__":
    main()
