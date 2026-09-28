"""Deploys busco_gca_flow to a slurm-cli work pool.

The worker clones --branch from GitHub at run time (same convention as
gb_prefect/deployments/deploy_flows.py), so commit and push the branch first --
including gb_prefect/flows/gb_busco_gca.py and gb_prefect/tasks/busco_gca.py.

Example (run from the repository root, with PREFECT_API_URL pointing at the
Prefect server):

    python gb_prefect/deployments/deploy_busco_gca.py \
        --branch dev/gb_prefect \
        --pool codon-slurm-smoke-pool \
        --outdir /path/to/output \
        --enscode /path/to/enscode \
        --asm-venv /path/to/asm_venv

This targets a slurm-cli pool (job_variables: cpu/memory/time_limit/partition/
setup_commands), unlike gb_prefect/deployments/deploy_flows.py's `process`-pool
flows -- kept separate since the parameter/job_variable shapes differ.
"""

import argparse

from prefect import flow  # type: ignore
from prefect.runner.storage import GitRepository  # type: ignore

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--branch", default="main", help="Git branch the worker clones.")
    parser.add_argument(
        "--pool",
        default="codon-slurm-smoke-pool",
        help="slurm-cli work pool to deploy to.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Base output directory (a <gca> subdirectory is created under it per run).",
    )
    parser.add_argument("--enscode", required=True, help="Path to the ENSCODE directory.")
    parser.add_argument(
        "--asm-venv",
        required=True,
        help="Path to the venv with gb_prefect + prefect installed; sourced by the worker "
        "before the flow run starts.",
    )
    parser.add_argument("--cpu", type=int, default=2)
    parser.add_argument("--memory", type=int, default=8, help="Memory in GB.")
    parser.add_argument("--time-limit", type=int, default=4, help="Wall time in hours.")
    parser.add_argument("--partition", default=None, help="Slurm partition, optional.")
    args = parser.parse_args()

    source = GitRepository(url="https://github.com/Ensembl/ensembl-genes-metadata.git", branch=args.branch)

    deployment_id = flow.from_source(
        source=source,
        entrypoint="gb_prefect/flows/gb_busco_gca.py:busco_gca_flow",
    ).deploy(
        name="busco-gca",
        work_pool_name=args.pool,
        parameters={"outdir": args.outdir, "enscode": args.enscode},
        job_variables={
            "cpu": args.cpu,
            "memory": args.memory,
            "time_limit": args.time_limit,
            "partition": args.partition,
            "working_dir": args.outdir,
            "setup_commands": [
                "module load nextflow/24.10.3",
                f"source {args.asm_venv}/bin/activate",
            ],
        },
        tags=["genebuild", "busco", args.branch],
        description=f"BUSCO for a single GCA from branch {args.branch}.",
    )
    print(f"Deployed busco-gca ({deployment_id})")
