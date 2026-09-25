"""Deploy Prefect flows

The worker clones --branch from GitHub at run time, so commit and push the branch first:
the flow code is taken from GitHub.

Credentials are not deployment parameters; the flows load them from Prefect Secret blocks
(create them once with gb_prefect/deployments/create_secrets.py).

Example (run from the repository root, with PREFECT_API_URL pointing at the Prefect server):

    python gb_prefect/deployments/deploy_flows.py \
        --branch main \
        --outdir /path/to/output \
        --work-pool <work_pool_name> \
        --enscode /path/to/enscode \
        --asm-venv /path/to/asm_venv \
        [--flows name_flow_1 name_flow2]

"""

import argparse

from prefect import flow  # type: ignore
from prefect.runner.storage import GitRepository  # type: ignore


# (entrypoint relative to the repo root, base deployment name, tag)
FLOWS = {
    "update": (
        "gb_prefect/flows/gb_registry_update.py:gb_registry_update_flow",
        "gb_metadata_update",
        "metadata_update",
    ),
    "by_date": (
        "gb_prefect/flows/gb_registry_update_by_date.py:gb_registry_update_by_date_flow",
        "gb_metadata_update_per_date",
        "metadata_update",
    ),
    "register": (
        "gb_prefect/flows/gb_registry.py:gb_registry_flow",
        "gb_metadata_registry",
        "metadata_registry",
    ),
    "register_by_gca": (
        "gb_prefect/flows/gb_registry_by_gca.py:gb_registry_by_gca_flow",
        "gb_metadata_registry_per_gca",
        "metadata_registry",
    ),
}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--branch", default="main", help="Git branch the worker clones.")
    parser.add_argument(
        "--name-suffix",
        default="",
        help="Appended to the deployment names; use '' to replace the production deployments.",
    )
    parser.add_argument("--flows", nargs="+", choices=sorted(FLOWS), default=sorted(FLOWS))
    parser.add_argument("--work-pool", default="genebuild-pool")
    parser.add_argument("--outdir", required=True, help="Base output directory.")
    parser.add_argument("--enscode", required=True, help="Path to the ENSCODE directory.")
    parser.add_argument("--asm-venv", required=True, help="Path to the virtual environment.")
    args = parser.parse_args()

    source = GitRepository(url="https://github.com/Ensembl/ensembl-genes-metadata.git", branch=args.branch)
    default_parameters = {"outdir": args.outdir, "enscode": args.enscode, "asm_venv": args.asm_venv}

    for key in args.flows:
        entrypoint, base_name, tag = FLOWS[key]
        deployment_id = flow.from_source(source=source, entrypoint=entrypoint).deploy(
            name=f"{base_name}{args.name_suffix}",
            work_pool_name=args.work_pool,
            parameters=default_parameters,
            tags=["genebuild", tag, args.branch],
            description=f"Assembly metadata {key} from branch {args.branch}.",
        )
        print(f"Deployed {entrypoint} as {base_name}{args.name_suffix} ({deployment_id})")
