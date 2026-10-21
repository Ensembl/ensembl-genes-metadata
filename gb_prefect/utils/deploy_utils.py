"""Deploy code for Nextflow flows on slurm-cli pools (used by gb_prefect/deployments/deploy_flows.py).

A deploy script lists its work pools (key -> default pool name) and its flows as FlowDeployment
entries, and calls main(); this module adds the command-line options -- including one
--<key>-pool per pool, to deploy to a different pool -- builds each deployment's parameters and
slurm-cli job variables, and deploys them (or prints them with --dry-run).

The flow code always comes from this repository (cloned by the worker at run time), whichever
repository the Nextflow pipeline itself lives in -- the pipeline is found under ENSCODE.
"""

import argparse
import json
import sys
from dataclasses import dataclass, field
from typing import Any, Dict, List, Mapping, Optional, Sequence

FLOW_REPO_URL = "https://github.com/Ensembl/ensembl-genes-metadata.git"
DEFAULT_NEXTFLOW_MODULE = "nextflow/25.04.6"


@dataclass(frozen=True)
class FlowDeployment:  # pylint: disable=too-many-instance-attributes
    """One deployment of a Nextflow flow to a slurm-cli work pool.

    Attributes:
        name: Deployment name (--name-suffix is appended).
        entrypoint: Flow entrypoint, relative to the repository root ("path.py:function").
        work_pool: Key of the flow's slurm-cli work pool in the script's pools; the pool name is
            that key's default, or --<key>-pool, or --work-pool (every flow).
        description: Shown in the Prefect UI; " from branch <branch>." is appended.
        tags: Extra tags, between "genebuild" and the branch tag.
        cpu, memory, time_limit: Slurm resources of the flow run's job -- the Nextflow head job,
            which must outlive the whole pipeline. memory in GB, time_limit in hours.
        concurrency_limit: Max simultaneous runs of this deployment (None: no limit).
        parameters: Extra default parameters, on top of outdir and enscode.
        separate_outdir: The flow's output normally lives apart from the shared --outdir (e.g.
            BUSCO's per-GCA run directories): a warning is printed if it isn't given its own
            with --flow-outdir.
    """

    name: str
    entrypoint: str
    work_pool: str
    description: str
    tags: Sequence[str] = ()
    cpu: int = 1
    memory: int = 4
    time_limit: int = 2
    concurrency_limit: Optional[int] = None
    parameters: Mapping[str, Any] = field(default_factory=dict)
    separate_outdir: bool = False


def _pool_dest(key: str) -> str:
    return f"pool_{key.replace('-', '_')}"


def parse_args(
    flows: Mapping[str, FlowDeployment], pools: Mapping[str, str], description: str
) -> argparse.Namespace:
    """Parse the deploy options; --flows chooses among flows' keys, --<key>-pool renames a pool.

    args.pools is set to the resolved pool names, by key.
    """
    unknown = sorted({flow.work_pool for flow in flows.values()} - set(pools))
    if unknown:
        raise ValueError(f"Flows use work pool key(s) missing from pools: {unknown}")

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--branch", default="main", help="Git branch the worker clones.")
    parser.add_argument(
        "--name-suffix",
        default="",
        help="Appended to the deployment names, e.g. _test; leave empty to replace the production ones.",
    )
    parser.add_argument("--flows", nargs="+", choices=sorted(flows), default=sorted(flows))
    parser.add_argument(
        "--outdir",
        required=True,
        help="Base output directory of every flow (also the job's working dir), unless --flow-outdir.",
    )
    parser.add_argument(
        "--flow-outdir",
        action="append",
        default=[],
        metavar="KEY=PATH",
        help="Base output directory of one flow instead of --outdir, e.g. busco_gca=/path; repeatable.",
    )
    parser.add_argument("--enscode", required=True, help="Path to the ENSCODE directory.")
    parser.add_argument(
        "--asm-venv",
        required=True,
        help="Venv with gb_prefect + prefect installed; activated by the worker before the flow run.",
    )
    parser.add_argument(
        "--nextflow-module",
        default=DEFAULT_NEXTFLOW_MODULE,
        help="Environment module loaded by the worker before the flow run.",
    )
    for key, default in pools.items():
        users = ", ".join(sorted(name for name, flow in flows.items() if flow.work_pool == key)) or "none"
        parser.add_argument(
            f"--{key}-pool",
            dest=_pool_dest(key),
            default=default,
            metavar="POOL",
            help=f"slurm-cli work pool of the {key} flows ({users}). Default: {default}.",
        )
    parser.add_argument(
        "--work-pool", help="slurm-cli work pool for every selected flow, e.g. a personal test pool."
    )
    parser.add_argument("--partition", help="Slurm partition, optional.")
    parser.add_argument("--cpu", type=int, help="Override every flow's CPU count.")
    parser.add_argument("--memory", type=int, help="Override every flow's memory, in GB.")
    parser.add_argument("--time-limit", type=int, help="Override every flow's wall time, in hours.")
    parser.add_argument(
        "--dry-run", action="store_true", help="Print the deployments that would be made; deploy nothing."
    )
    args = parser.parse_args()

    flow_outdirs: Dict[str, str] = {}
    for value in args.flow_outdir:
        key, sep, path = value.partition("=")
        if not sep or not path:
            parser.error(f"--flow-outdir must be KEY=PATH, got {value!r}")
        if key not in flows:
            parser.error(f"--flow-outdir: unknown flow {key!r} (choose from {', '.join(sorted(flows))})")
        flow_outdirs[key] = path
    args.flow_outdir = flow_outdirs
    args.pools = {key: getattr(args, _pool_dest(key)) for key in pools}
    return args


def build_deployment(
    flow_deployment: FlowDeployment, args: argparse.Namespace, outdir: Optional[str] = None
) -> Dict[str, Any]:
    """Return the deploy() arguments of flow_deployment for the given options.

    outdir is the flow's base output directory (its --flow-outdir); defaults to args.outdir.
    """
    outdir = outdir or args.outdir

    def override(value: Optional[int], default: int) -> int:
        return default if value is None else value

    spec: Dict[str, Any] = {
        "name": f"{flow_deployment.name}{args.name_suffix}",
        "entrypoint": flow_deployment.entrypoint,
        "work_pool_name": args.work_pool or args.pools[flow_deployment.work_pool],
        "parameters": {"outdir": outdir, "enscode": args.enscode, **flow_deployment.parameters},
        "job_variables": {
            "cpu": override(args.cpu, flow_deployment.cpu),
            "memory": override(args.memory, flow_deployment.memory),
            "time_limit": override(args.time_limit, flow_deployment.time_limit),
            "partition": args.partition,
            "working_dir": outdir,
            "setup_commands": [
                f"module load {args.nextflow_module}",
                f"source {args.asm_venv}/bin/activate",
            ],
        },
        "tags": ["genebuild", *flow_deployment.tags, args.branch],
        "description": f"{flow_deployment.description} from branch {args.branch}.",
    }
    if flow_deployment.concurrency_limit is not None:
        spec["concurrency_limit"] = flow_deployment.concurrency_limit
    return spec


def main(flows: Mapping[str, FlowDeployment], pools: Mapping[str, str], description: str) -> None:
    """Parse the options, then deploy (or with --dry-run print) the selected flows."""
    args = parse_args(flows, pools, description)
    for key in args.flows:
        if flows[key].separate_outdir and key not in args.flow_outdir:
            print(
                f"WARNING: {key} uses --outdir {args.outdir}; its output usually has its own "
                f"directory (--flow-outdir {key}=PATH).",
                file=sys.stderr,
            )
    specs: List[Dict[str, Any]] = [
        build_deployment(flows[key], args, args.flow_outdir.get(key)) for key in args.flows
    ]

    if args.dry_run:
        print(json.dumps(specs, indent=2))
        print(f"Dry run: {len(specs)} deployment(s) resolved, nothing deployed.")
        return

    # Imported here so --dry-run works without a Prefect server configured.
    from prefect import flow  # type: ignore  # pylint: disable=import-outside-toplevel
    from prefect.runner.storage import GitRepository  # type: ignore  # pylint: disable=import-outside-toplevel

    source = GitRepository(url=FLOW_REPO_URL, branch=args.branch)
    for spec in specs:
        entrypoint = spec.pop("entrypoint")
        # from_source is typed sync-or-async; called from sync code it returns the Flow.
        remote_flow = flow.from_source(source=source, entrypoint=entrypoint)
        deployment_id = remote_flow.deploy(**spec)  # type: ignore[union-attr]
        print(f"Deployed {entrypoint} as {spec['name']} on {spec['work_pool_name']} ({deployment_id})")
