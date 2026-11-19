from datetime import datetime
from pathlib import Path
from typing import Optional

from prefect import flow  # type: ignore

from gb_prefect.models.pipeline_options import TaskRunOptions
from gb_prefect.pipelines.assembly_metadata_update import asm_update_run
from gb_prefect.tasks.nextflow import run_nextflow_pipeline
from gb_prefect.tasks.registry_query import DEFAULT_GCA_QUERY, collect_gca_list
from gb_prefect.utils.credentials_utils import PipelineCredentials, resolve_credentials
from gb_prefect.utils.enscode_utils import resolve_enscode
from gb_prefect.utils.path_utils import run_folder

GCA_LIST_FILE = "gca_list.txt"


@flow(name="gb_registry_update_from_query", log_prints=True)
def gb_registry_update_from_query_flow(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    outdir: str,
    query: str = DEFAULT_GCA_QUERY,
    slack_report: bool = True,
    folder_name: Optional[str] = None,
    enscode: Optional[str] = None,
    dry_run: bool = False,
    credentials: Optional[PipelineCredentials] = None,
    metadata_secret_block: Optional[str] = None,
    slack_secret_block: Optional[str] = None,
):
    """Collect a GCA list from the metadata DB with a SQL query, then run the assembly
    metadata update Nextflow pipeline on it (as gb_registry_update_flow does for a file).

    query: one SELECT statement whose first column is the GCA accession. Default: every GCA
        in genebuild_status not handed_over, live or coming_soon (DEFAULT_GCA_QUERY). It runs
        in a read-only transaction.

    The list is written to {outdir}/{folder_name}/gca_list.txt, one accession per line;
    folder_name defaults to asm_update_query_<today>. If the query returns no GCAs the
    pipeline isn't run.

    dry_run: the query still runs (it only reads, and the list is worth checking), so the
    metadata credentials are loaded; only the pipeline isn't executed.
    """
    # Resolve credentials: the query needs the real DB params, even on a dry run.
    credentials = resolve_credentials(
        credentials,
        slack_report=slack_report,
        metadata_secret_block=metadata_secret_block,
        slack_secret_block=slack_secret_block,
    )
    run_dir = run_folder(outdir, folder_name, f"asm_update_query_{datetime.now().strftime('%Y-%m-%d')}")
    gca_list = str(Path(run_dir) / GCA_LIST_FILE)

    # Collect the GCAs from the registry.
    gcas = collect_gca_list(query, credentials.metadata_db_params, gca_list)
    if not gcas:
        print("The query returned no GCAs; nothing to update.")
        return None

    # Run the assembly metadata update pipeline on the collected list.
    run = asm_update_run(
        outdir=run_dir,
        enscode=resolve_enscode(enscode, dry_run),
        metadata_params_string=credentials.metadata_params_string,
        gca_list=gca_list,
        slack_report=slack_report,
        slack_params=credentials.slack_params,
    )
    return run_nextflow_pipeline(run, TaskRunOptions(dry_run=dry_run))
