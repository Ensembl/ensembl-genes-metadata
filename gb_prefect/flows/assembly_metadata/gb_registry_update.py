from datetime import datetime
from typing import Optional

from prefect import flow  # type: ignore

from gb_prefect.models.pipeline_options import TaskRunOptions
from gb_prefect.pipelines.assembly_metadata_update import asm_update_run
from gb_prefect.tasks.nextflow import run_nextflow_pipeline
from gb_prefect.utils.credentials_utils import PipelineCredentials, resolve_credentials
from gb_prefect.utils.enscode_utils import resolve_enscode
from gb_prefect.utils.path_utils import check_input_file, run_folder


@flow(name="gb_registry_update", log_prints=True)
def gb_registry_update_flow(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    gca_list: str,
    outdir: str,
    slack_report: bool = True,
    folder_name: Optional[str] = None,
    enscode: Optional[str] = None,
    dry_run: bool = False,
    credentials: Optional[PipelineCredentials] = None,
    metadata_secret_block: Optional[str] = None,
    slack_secret_block: Optional[str] = None,
):
    """Run the assembly metadata update Nextflow pipeline for the GCA list (the whole file, as
    given). See gb_registry_update_by_date_flow for screening by release date instead.

    gca_list must be an absolute path. Output goes to {outdir}/{folder_name}; folder_name
    defaults to asm_update_<today>.
    """
    # Resolve credentials
    credentials = resolve_credentials(
        credentials,
        slack_report=slack_report,
        dry_run=dry_run,
        metadata_secret_block=metadata_secret_block,
        slack_secret_block=slack_secret_block,
    )

    # Run the assembly metadata update pipeline for the specified GCA accessions.
    run = asm_update_run(
        outdir=run_folder(outdir, folder_name, f"asm_update_{datetime.now().strftime('%Y-%m-%d')}"),
        enscode=resolve_enscode(enscode, dry_run),
        metadata_params_string=credentials.metadata_params_string,
        gca_list=check_input_file(gca_list, dry_run),
        slack_report=slack_report,
        slack_params=credentials.slack_params,
    )
    return run_nextflow_pipeline(run, TaskRunOptions(dry_run=dry_run))
