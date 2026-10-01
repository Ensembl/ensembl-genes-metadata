"""Deploys busco_dispatch_flow (the BUSCO orchestrator, scenario 1: explicit CSV input).

The worker clones --branch from GitHub at run time (same convention as
gb_prefect/deployments/deploy_flows.py), so commit and push the branch first.

Credentials are not a deployment parameter -- the flow loads them from the
gb-metadata-db-params Prefect Secret block (create it once with
gb_prefect/deployments/create_secrets.py).

csv_file is left unset here -- pass it per run, e.g.:

    prefect deployment run 'BUSCO_dispatch/busco-dispatch' --param csv-file=/path/to/gcas.csv

Example (run from the repository root, with PREFECT_API_URL pointing at the Prefect server):

    python gb_prefect/deployments/deploy_busco_dispatch.py --branch dev/gb_prefect
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
    args = parser.parse_args()

    source = GitRepository(url="https://github.com/Ensembl/ensembl-genes-metadata.git", branch=args.branch)

    deployment_id = flow.from_source(
        source=source,
        entrypoint="gb_prefect/flows/gb_busco_dispatch.py:busco_dispatch_flow",
    ).deploy(
        name=f"busco-dispatch{args.name_suffix}",
        work_pool_name=args.work_pool,
        tags=["genebuild", "busco", args.branch],
        description=f"BUSCO orchestrator (explicit CSV input) from branch {args.branch}.",
    )
    print(f"Deployed busco-dispatch ({deployment_id})")
