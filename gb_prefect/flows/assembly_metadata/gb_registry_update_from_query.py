import hashlib
from pathlib import Path
from typing import Any, Dict, List, Optional, cast

from prefect import flow, task  # type: ignore
from prefect.client.schemas.objects import FlowRun  # type: ignore
from prefect.deployments import run_deployment  # type: ignore
from prefect.runtime import flow_run  # type: ignore

from gb_prefect.tasks.registry_query import (
    DEFAULT_GCA_QUERY,
    collect_gca_list,
    read_gca_file,
    split_batches,
    write_gca_file,
)
from gb_prefect.utils.credentials_utils import PipelineCredentials, resolve_credentials
from gb_prefect.utils.path_utils import run_folder
from gb_prefect.utils.run_utils import is_retry, run_name_suffix, run_start

GCA_LIST_FILE = "gca_list.txt"
BATCH_SIZE = 1000
UPDATE_DEPLOYMENT = "gb_registry_update/gb_metadata_update"
# How often to check whether a batch's run has finished.
BATCH_POLL_SECONDS = 60


def batch_idempotency_key(batch_name: str, gcas: List[str]) -> str:
    """Key for triggering one batch: the same flow run, batch and GCAs give the same key, so a
    retry of this flow gets the batch's existing run back instead of starting it again."""
    digest = hashlib.sha1("\n".join(gcas).encode("utf-8")).hexdigest()[:12]
    return f"{flow_run.id}-{batch_name}-{digest}"


@task(log_prints=True)
def run_update_batch(
    deployment_name: str, parameters: Dict[str, Any], idempotency_key: str
) -> Dict[str, Any]:
    """Trigger one batch's update run and wait until it finishes. Never raises: a failure to
    trigger is returned as state NOT_STARTED, so the next batches still run."""
    try:
        # run_deployment is typed sync-or-async; called from sync code it returns the FlowRun.
        child = cast(
            FlowRun,
            run_deployment(
                name=deployment_name,
                parameters=parameters,
                idempotency_key=idempotency_key,
                poll_interval=BATCH_POLL_SECONDS,
            ),
        )
    except Exception as err:  # pylint: disable=broad-exception-caught
        print(f"Could not trigger {deployment_name}: {err}")
        return {"flow_run_id": None, "state": "NOT_STARTED", "error": str(err)}
    state = child.state.type.value if child.state else "UNKNOWN"
    print(f"Run {child.id} ({child.name}) finished: {state}")
    return {"flow_run_id": str(child.id), "state": state}


@flow(name="gb_registry_update_from_query", log_prints=True)
def gb_registry_update_from_query_flow(  # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals
    outdir: str,
    query: str = DEFAULT_GCA_QUERY,
    batch_size: int = BATCH_SIZE,
    slack_report: bool = True,
    folder_name: Optional[str] = None,
    enscode: Optional[str] = None,
    dry_run: bool = False,
    deployment_name: str = UPDATE_DEPLOYMENT,
    credentials: Optional[PipelineCredentials] = None,
    metadata_secret_block: Optional[str] = None,
    slack_secret_block: Optional[str] = None,
) -> List[Dict[str, Any]]:
    """Collect a GCA list from the metadata DB with a SQL query, split it into batches of
    batch_size, and run the assembly metadata update on each batch, one after another: each
    batch is its own run of deployment_name (gb_registry_update_flow), and the next batch is
    triggered only once the previous run has finished.

    query: one SELECT statement whose first column is the GCA accession. Default: every GCA
        in genebuild_status not handed_over, live or coming_soon (DEFAULT_GCA_QUERY). It runs
        in a read-only transaction.

    Output, under {outdir}/{folder_name} (folder_name defaults to
    asm_update_query_<run start date>_<flow run name>, which a retry keeps):
        gca_list.txt                the whole list, one accession per line
        batch_001/gca_list.txt      each batch's list; the batch's run writes its pipeline
        batch_002/...               output in the same folder

    A batch whose run fails (or can't be triggered) doesn't stop the others; once every batch
    has run, this flow fails if any of them did, naming them. Retry a failed batch from its own
    run.

    Retrying this flow doesn't run the query again: the list in gca_list.txt is reused (the
    query runs only if the first attempt didn't get as far as writing it), so the batches are
    the same, and each batch already triggered gets its existing run back instead of a new
    one. A new run of the flow (not a retry) always runs the query, even into an existing
    folder_name.

    Each batch's run sends its own Slack report when slack_report is true. enscode,
    credentials and the Secret block names are passed on to each run when given; otherwise
    the update deployment's own defaults apply.

    dry_run: the query still runs (it only reads) and the list and batch files are written,
    but no batch run is triggered.
    """
    # The default folder comes from the flow run, so a retry (even on a later day) finds it again.
    default_folder = f"asm_update_query_{run_start().strftime('%Y-%m-%d')}{run_name_suffix()}"
    run_dir = Path(run_folder(outdir, folder_name, default_folder))
    gca_list = run_dir / GCA_LIST_FILE

    # Collect the GCAs from the registry -- or, on a retry, reuse the list the first attempt
    # collected, so the batches don't change under runs already triggered.
    if is_retry() and gca_list.is_file():
        gcas = read_gca_file(str(gca_list))
        print(f"Retry: reusing the {len(gcas)} GCA(s) in {gca_list}; the query is not run again.")
    else:
        # The query needs the real DB params, even on a dry run.
        db_params = resolve_credentials(
            credentials, metadata_secret_block=metadata_secret_block
        ).metadata_db_params
        gcas = collect_gca_list(query, db_params, str(gca_list))
    if not gcas:
        print("The query returned no GCAs; nothing to update.")
        return []
    batches = split_batches(gcas, batch_size)
    print(f"{len(gcas)} GCA(s) in {len(batches)} batch(es) of up to {batch_size}.")

    # Run the batches one after another.
    results = []
    for number, batch in enumerate(batches, start=1):
        batch_name = f"batch_{number:03d}"
        parameters: Dict[str, Any] = {
            "gca_list": write_gca_file(str(run_dir / batch_name / GCA_LIST_FILE), batch),
            "outdir": str(run_dir),
            "folder_name": batch_name,
            "slack_report": slack_report,
        }
        optional = {
            "enscode": enscode,
            "credentials": credentials,
            "metadata_secret_block": metadata_secret_block,
            "slack_secret_block": slack_secret_block,
        }
        parameters.update({key: value for key, value in optional.items() if value is not None})
        print(f"{batch_name} ({number}/{len(batches)}): {len(batch)} GCA(s), {batch[0]} .. {batch[-1]}")

        if dry_run:
            print(f"Dry run: not triggering {deployment_name} for {batch_name}.")
            result: Dict[str, Any] = {"flow_run_id": None, "state": "DRY_RUN"}
        else:
            result = run_update_batch(
                deployment_name, parameters, batch_idempotency_key(batch_name, batch)
            )
        results.append({"batch": batch_name, "gca_count": len(batch), **result})

    # Report every batch, then fail if any of them didn't complete.
    for result in results:
        print(
            f"{result['batch']}: {result['gca_count']} GCA(s), {result['state']}, "
            f"run {result['flow_run_id']}"
        )
    failed = [r["batch"] for r in results if r["state"] not in ("COMPLETED", "DRY_RUN")]
    if failed:
        raise RuntimeError(f"{len(failed)} of {len(results)} batch(es) didn't complete: {', '.join(failed)}")
    return results
