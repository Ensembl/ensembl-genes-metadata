"""Shared code of the deploy scripts in gb_prefect/deployments/ for Nextflow flows on slurm-cli pools.

A deploy script lists its flows as FlowDeployment entries and calls main(); this module adds the
common command-line options, builds each deployment's parameters and slurm-cli job variables,
and deploys them (or prints them with --dry-run).

The flow code always comes from this repository (cloned by the worker at run time), whichever
repository the Nextflow pipeline itself lives in -- the pipeline is found under ENSCODE.
"""

import argparse
import json
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
        work_pool: Default slurm-cli work pool (--work-pool overrides it for every flow).
        description: Shown in the Prefect UI; " from branch <branch>." is appended.
        tags: Extra tags, between "genebuild" and the branch tag.
        cpu, memory, time_limit: Slurm resources of the flow run's job -- the Nextflow head job,
            which must outlive the whole pipeline. memory in GB, time_limit in hours.
        concurrency_limit: Max simultaneous runs of this deployment (None: no limit).
        parameters: Extra default parameters, on top of outdir and enscode.
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


def parse_args(flows: Mapping[str, FlowDeployment], description: str) -> argparse.Namespace:
    """Parse the deploy options shared by every deploy script; --flows chooses among flows' keys."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--branch", default="main", help="Git branch the worker clones.")
    parser.add_argument(
        "--name-suffix",
        default="",
        help="Appended to the deployment names, e.g. _test; leave empty to replace the production ones.",
    )
    parser.add_argument("--flows", nargs="+", choices=sorted(flows), default=sorted(flows))
    parser.add_argument("--outdir", required=True, help="Base output directory (also the job's working dir).")
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
    parser.add_argument("--work-pool", help="slurm-cli work pool for every flow (default: each flow's own).")
    parser.add_argument("--partition", help="Slurm partition, optional.")
    parser.add_argument("--cpu", type=int, help="Override every flow's CPU count.")
    parser.add_argument("--memory", type=int, help="Override every flow's memory, in GB.")
    parser.add_argument("--time-limit", type=int, help="Override every flow's wall time, in hours.")
    parser.add_argument(
        "--dry-run", action="store_true", help="Print the deployments that would be made; deploy nothing."
    )
    return parser.parse_args()


def build_deployment(flow_deployment: FlowDeployment, args: argparse.Namespace) -> Dict[str, Any]:
    """Return the deploy() arguments of flow_deployment for the given options."""

    def override(value: Optional[int], default: int) -> int:
        return default if value is None else value

    spec: Dict[str, Any] = {
        "name": f"{flow_deployment.name}{args.name_suffix}",
        "entrypoint": flow_deployment.entrypoint,
        "work_pool_name": args.work_pool or flow_deployment.work_pool,
        "parameters": {"outdir": args.outdir, "enscode": args.enscode, **flow_deployment.parameters},
        "job_variables": {
            "cpu": override(args.cpu, flow_deployment.cpu),
            "memory": override(args.memory, flow_deployment.memory),
            "time_limit": override(args.time_limit, flow_deployment.time_limit),
            "partition": args.partition,
            "working_dir": args.outdir,
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


def main(flows: Mapping[str, FlowDeployment], description: str) -> None:
    """Parse the options, then deploy (or with --dry-run print) the selected flows."""
    args = parse_args(flows, description)
    specs: List[Dict[str, Any]] = [build_deployment(flows[key], args) for key in args.flows]

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
