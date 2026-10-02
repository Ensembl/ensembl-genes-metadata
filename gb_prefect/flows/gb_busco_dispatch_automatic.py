import argparse
import json
import logging
from typing import Optional

from prefect import flow  # type: ignore

from gb_metadata.update_busco_events import (
    get_pending_candidates,
    insert_busco_candidates,
    mark_existing_busco_done,
)
from gb_prefect.models.pipeline_options import PipelineCredentials
from gb_prefect.tasks.busco_dispatch import dispatch_rows, get_free_slots
from gb_prefect.utils.credentials_utils import (
    DEFAULT_METADATA_SECRET_BLOCK,
    DEFAULT_SLACK_SECRET_BLOCK,
    resolve_credentials,
)

BUSCO_GCA_DEPLOYMENT = "BUSCO_gca/busco-gca"
BUSCO_POOL = "codon-slurm-smoke-pool"


@flow(name="BUSCO_dispatch_automatic", log_prints=True)
def busco_dispatch_automatic_flow(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    pool_name: str = BUSCO_POOL,
    dry_run: bool = False,
    deployment_name: str = BUSCO_GCA_DEPLOYMENT,
    credentials: Optional[PipelineCredentials] = None,
    metadata_secret_block: str = DEFAULT_METADATA_SECRET_BLOCK,
):
    """Automatically select and trigger busco_gca_flow for pending candidates, sized to
    pool_name's currently free concurrency slots.

    Registry update is always on here (no option to turn it off, unlike scenarios 1/2):
    every picked GCA has genome_busco.status set to in_progress right after a successful
    trigger, so the next automatic run doesn't pick it again while it's still running, and
    each triggered run gets the Nextflow pipeline's --update_registry, which marks it done.

    Each triggered run also gets post_run_cleanup=True (scenario 3 only): a failed run sets
    genome_busco.status to failed and keeps only its logs; a successful run's directory is
    removed. See busco_gca_flow.

    Refreshes assembly_events first (same as scenarios 1/2). Selects up to N pending
    candidates (genome_busco.status=pending), ordered high -> medium -> low priority;
    large_genome is excluded -- that tier is run manually, never picked here. N is the
    pool's free slots (concurrency_limit - active_slots); refuses to run at all if the pool
    has no concurrency_limit set, rather than guessing a fallback cap.

    No --force here: automatic mode only ever picks genuinely pending candidates, never
    re-triggers something already done.

    Manual trigger only for now -- not yet wired to a schedule or to a "pool has space"
    event trigger.
    """
    file_handler = logging.FileHandler("busco_dispatch_automatic.log", mode="w")
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

    free_slots = get_free_slots(pool_name)
    if free_slots == 0:
        print(f"No free slots in '{pool_name}'; nothing to dispatch.")
        return {"triggered": [], "skipped": []}

    candidates = get_pending_candidates(db_params, free_slots)
    print(f"Selected {len(candidates)} candidate(s) for {free_slots} free slot(s)")

    rows = [
        {"gca": gca, "taxon_id": str(taxon_id), "busco_dataset": ""}
        for gca, taxon_id, _priority in candidates
    ]

    return dispatch_rows(
        rows,
        db_params,
        force=False,
        dry_run=dry_run,
        deployment_name=deployment_name,
        update_registry=True,
        post_run_cleanup=True,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pool-name", default=BUSCO_POOL, help="slurm-cli work pool to size the batch against."
    )
    parser.add_argument(
        "--metadata-params-string",
        required=False,
        help="JSON string with metadata DB connection parameters. If omitted, credentials are "
        "loaded from the Prefect Secret blocks (see gb_prefect/deployments/create_secrets.py).",
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

    busco_dispatch_automatic_flow(
        pool_name=args.pool_name,
        dry_run=args.dry_run,
        deployment_name=args.deployment_name,
        credentials=(
            PipelineCredentials(metadata_params_string=args.metadata_params_string)
            if args.metadata_params_string
            else None
        ),
    )
