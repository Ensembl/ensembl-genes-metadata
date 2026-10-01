"""Deploys busco_dispatch_flow (the BUSCO orchestrator, scenario 1: explicit CSV input).

The worker clones --branch from GitHub at run time (same convention as
gb_prefect/deployments/deploy_flows.py), so commit and push the branch first.

Credentials are not a deployment parameter -- the flow loads them from the
gb-metadata-db-params Prefect Secret block (create it once with
gb_prefect/deployments/create_secrets.py).

csv_file is left unset here -- pass it per run, e.g.:

    prefect deployment run 'BUSCO_dispatch/busco-dispatch' --param csv-file=/path/to/gcas.csv

A `process`-type pool has no setup_commands/venv-activation mechanism like slurm-cli does --
the flow run just inherits whatever Python environment the worker process polling that pool
happens to be running in. If that's not the venv you want (e.g. you need your own, up-to-date
checkout rather than whatever a shared worker is using), pass --asm-venv: it overrides the
job's `command` job_variable to invoke that venv's python directly, instead of the worker's own
sys.executable (Prefect only substitutes its own default if `command` is left unset).

Example (run from the repository root, with PREFECT_API_URL pointing at the Prefect server):

    python gb_prefect/deployments/deploy_busco_dispatch.py --branch dev/gb_prefect \
        --work-pool genebuild-pool --asm-venv /path/to/your/venv
"""

import argparse

from prefect import flow  # type: ignore
from prefect.runner.storage import GitRepository  # type: ignore

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--branch",
        required=True,
        help="Git branch the worker clones."
    )
    parser.add_argument(
        "--name-suffix",
        default="",
        help="Appended to the deployment name; use '' to replace the production deployment.",
    )
    parser.add_argument(
        "--work-pool",
        required=True,
        help="process-type work pool to deploy to."
    )
    parser.add_argument(
        "--asm-venv",
        default=None,
        help="Path to a venv with gb_prefect + gb_metadata + prefect installed, e.g. your own "
        "rather than whatever venv the pool's worker process happens to be running in. "
        "Overrides the job's command to use that venv's python directly.",
    )
    args = parser.parse_args()

    source = GitRepository(url="https://github.com/Ensembl/ensembl-genes-metadata.git", branch=args.branch)

    job_variables = {}
    if args.asm_venv:
        job_variables["command"] = f"{args.asm_venv.rstrip('/')}/bin/python -m prefect.engine"

    deployment_id = flow.from_source(
        source=source,
        entrypoint="gb_prefect/flows/gb_busco_dispatch.py:busco_dispatch_flow",
    ).deploy(
        name=f"busco-dispatch{args.name_suffix}",
        work_pool_name=args.work_pool,
        job_variables=job_variables,
        tags=["genebuild", "busco", args.branch],
        description=f"BUSCO orchestrator (explicit CSV input) from branch {args.branch}.",
    )
    print(f"Deployed busco-dispatch ({deployment_id})")
