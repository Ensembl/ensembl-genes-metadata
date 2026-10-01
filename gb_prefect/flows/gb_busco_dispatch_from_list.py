import argparse
import json
import logging
from typing import List, Optional

from prefect import flow  # type: ignore

from gb_metadata.update_busco_events import get_taxon_id, insert_busco_candidates, mark_existing_busco_done
from gb_prefect.models.pipeline_options import PipelineCredentials
from gb_prefect.tasks.busco_dispatch import dispatch_rows
from gb_prefect.utils.credentials_utils import (
    DEFAULT_METADATA_SECRET_BLOCK,
    DEFAULT_SLACK_SECRET_BLOCK,
    resolve_credentials,
)

BUSCO_GCA_DEPLOYMENT = "BUSCO_gca/busco-gca"


def _read_gca_list_file(path: str) -> List[str]:
    """Read one GCA accession per line from a plain text file."""
    with open(path, encoding="utf-8") as f:
        return [line.strip() for line in f if line.strip()]


@flow(name="BUSCO_dispatch_from_list", log_prints=True)
def busco_dispatch_from_list_flow(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    gca_list: Optional[List[str]] = None,
    gca_file: Optional[str] = None,
    force: bool = False,
    dry_run: bool = False,
    deployment_name: str = BUSCO_GCA_DEPLOYMENT,
    credentials: Optional[PipelineCredentials] = None,
    metadata_secret_block: str = DEFAULT_METADATA_SECRET_BLOCK,
):
    """Trigger busco_gca_flow for each GCA given directly -- gca_list and/or gca_file (one
    accession per line); values from both are combined and de-duplicated.

    taxon_id is resolved per GCA from the metadata DB (assembly.lowest_taxon_id). A GCA not
    found in the DB at all is skipped -- not run -- since there's no taxon_id to dispatch
    with. Otherwise behaves exactly like busco_dispatch_flow (scenario 1): assembly_events is
    refreshed first, and a GCA already marked done is skipped unless force=True. The shared
    trigger logic lives in gb_prefect.tasks.busco_dispatch.dispatch_rows.

    credentials is optional: when omitted (the normal case for a deployment trigger), DB
    connection params are loaded from the metadata_secret_block Prefect Secret block instead.
    """
    file_handler = logging.FileHandler("busco_dispatch_from_list.log", mode="w")
    file_handler.setFormatter(logging.Formatter("%(asctime)s:%(levelname)s:%(message)s"))
    root_logger = logging.getLogger()
    root_logger.addHandler(file_handler)
    root_logger.setLevel(logging.INFO)

    credentials = resolve_credentials(
        credentials, metadata_secret_block, DEFAULT_SLACK_SECRET_BLOCK, slack_report=False
    )
    db_params = json.loads(credentials.metadata_params_string)

    print("Refreshing assembly_events before dispatching...")
    mark_existing_busco_done(db_params, execute=not dry_run)
    insert_busco_candidates(db_params, execute=not dry_run)

    gcas = list(gca_list or [])
    if gca_file:
        gcas += _read_gca_list_file(gca_file)
    seen = set()
    unique_gcas = [g for g in gcas if not (g in seen or seen.add(g))]
    print(f"Resolving taxon_id for {len(unique_gcas)} GCA(s)")

    rows = []
    not_found = []
    for gca in unique_gcas:
        taxon_id = get_taxon_id(gca, db_params)
        if taxon_id is None:
            print(f"Skipping {gca}: not found in the metadata DB")
            not_found.append(gca)
            continue
        rows.append({"gca": gca, "taxon_id": str(taxon_id), "busco_dataset": ""})

    result = dispatch_rows(rows, db_params, force, dry_run, deployment_name)
    result["not_found"] = not_found
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--gca", nargs="+", default=None, help="One or more GCA accessions.")
    parser.add_argument("--gca-file", default=None, help="Text file with one GCA accession per line.")
    parser.add_argument(
        "--metadata-params-string",
        required=False,
        help="JSON string with metadata DB connection parameters. If omitted, credentials are "
        "loaded from the Prefect Secret blocks (see gb_prefect/deployments/create_secrets.py).",
    )
    parser.add_argument(
        "--force", action="store_true", help="Trigger a GCA even if already marked done."
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Refresh assembly_events in log-only mode and don't trigger any runs.",
    )
    parser.add_argument(
        "--deployment-name",
        default=BUSCO_GCA_DEPLOYMENT,
        help="Deployment to trigger per GCA (flow_name/deployment_name).",
    )
    args = parser.parse_args()

    if not args.gca and not args.gca_file:
        parser.error("Pass at least one of --gca / --gca-file.")

    busco_dispatch_from_list_flow(
        gca_list=args.gca,
        gca_file=args.gca_file,
        force=args.force,
        dry_run=args.dry_run,
        deployment_name=args.deployment_name,
        credentials=(
            PipelineCredentials(metadata_params_string=args.metadata_params_string)
            if args.metadata_params_string
            else None
        ),
    )
