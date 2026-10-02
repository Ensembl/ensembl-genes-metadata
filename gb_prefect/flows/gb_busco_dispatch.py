import argparse
import json
import logging
from typing import Optional

from prefect import flow  # type: ignore

from gb_metadata.update_busco_events import insert_busco_candidates, mark_existing_busco_done
from gb_prefect.models.pipeline_options import PipelineCredentials
from gb_prefect.tasks.busco_dispatch import dispatch_rows
from gb_prefect.utils.credentials_utils import (
    DEFAULT_METADATA_SECRET_BLOCK,
    DEFAULT_SLACK_SECRET_BLOCK,
    resolve_credentials,
)
from gb_prefect.utils.stats_split_csv import read_gca_csv

BUSCO_GCA_DEPLOYMENT = "BUSCO_gca/busco-gca"


@flow(name="BUSCO_dispatch", log_prints=True)
def busco_dispatch_flow(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    csv_file: str,
    force: bool = False,
    dry_run: bool = False,
    update_registry: bool = False,
    deployment_name: str = BUSCO_GCA_DEPLOYMENT,
    credentials: Optional[PipelineCredentials] = None,
    metadata_secret_block: str = DEFAULT_METADATA_SECRET_BLOCK,
):
    """Trigger busco_gca_flow for each GCA listed in csv_file (gca,taxon_id,busco_dataset).

    Always refreshes assembly_events first (gb_metadata.update_busco_events), so status
    checks below are current. A GCA already marked genome_busco.status=done is skipped
    unless force=True -- every other row runs regardless of its status/priority, since it
    was asked for explicitly in the input CSV rather than selected automatically.

    update_registry=False (the default) just triggers; True also marks genome_busco.status
    as in_progress right after a successful trigger, and passes update_registry on to each
    triggered run so the Nextflow pipeline's --update_registry loads the results and marks
    it done on completion (gb_prefect.tasks.busco_dispatch.dispatch_rows).

    credentials is optional: when omitted (the normal case for a deployment trigger, where
    it's never set as a parameter), DB connection params are loaded from the
    metadata_secret_block Prefect Secret block instead (same resolve_credentials() helper
    gb_registry_update_by_date_flow uses), so credentials never appear as a visible
    deployment parameter. Pass credentials explicitly to override, e.g. for local standalone
    testing without a Prefect server/Secret block available.

    This is scenario 1 of the orchestrator (explicit CSV input). The shared trigger logic
    (status check, skip/force, run_deployment) lives in gb_prefect.tasks.busco_dispatch.
    dispatch_rows, also used by scenario 2 (gb_busco_dispatch_from_list.py, bare GCA list)
    and scenario 3 (gb_busco_dispatch_automatic.py, automatic candidate selection sized to
    free work pool slots, where update_registry is always on).
    """
    # logging.basicConfig() is a no-op here: importing `prefect` already attaches a
    # PrefectConsoleHandler to the root logger, so basicConfig's "only if no handlers
    # exist yet" guard silently skips. Attach a FileHandler explicitly instead, so
    # update_busco_events's logging.info() calls (queries, candidate counts) land in a
    # real file without removing Prefect's own console handler.
    file_handler = logging.FileHandler("busco_dispatch.log", mode="w")
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

    rows = read_gca_csv(csv_file)
    print(f"Read {len(rows)} GCA(s) from {csv_file}")

    return dispatch_rows(rows, db_params, force, dry_run, deployment_name, update_registry)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--csv-file", required=True, help="CSV with gca,taxon_id,busco_dataset columns."
    )
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
        "--update_registry",
        action="store_true",
        help="Mark genome_busco.status as in_progress right after a successful trigger.",
    )
    parser.add_argument(
        "--deployment-name",
        default=BUSCO_GCA_DEPLOYMENT,
        help="Deployment to trigger per GCA (flow_name/deployment_name).",
    )
    args = parser.parse_args()

    busco_dispatch_flow(
        csv_file=args.csv_file,
        force=args.force,
        dry_run=args.dry_run,
        update_registry=args.update_registry,
        deployment_name=args.deployment_name,
        credentials=(
            PipelineCredentials(metadata_params_string=args.metadata_params_string)
            if args.metadata_params_string
            else None
        ),
    )
