from datetime import datetime
from typing import Optional

from prefect import flow  # type: ignore

from gb_prefect.models.pipeline_options import PipelineCredentials, TaskRunOptions
from gb_prefect.pipelines.assembly_metadata_update import asm_update_run
from gb_prefect.tasks.nextflow import run_nextflow_pipeline
from gb_prefect.utils.credentials_utils import (
    DEFAULT_METADATA_SECRET_BLOCK,
    DEFAULT_SLACK_SECRET_BLOCK,
    resolve_credentials,
)
from gb_prefect.utils.enscode_utils import resolve_enscode
from gb_prefect.utils.path_utils import run_folder


@flow(name="gb_registry_update_by_date", log_prints=True)
def gb_registry_update_by_date_flow(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    screen_date: str,
    outdir: str,
    slack_report: bool = True,
    folder_name: Optional[str] = None,
    enscode: Optional[str] = None,
    dry_run: bool = False,
    credentials: Optional[PipelineCredentials] = None,
    metadata_secret_block: str = DEFAULT_METADATA_SECRET_BLOCK,
    slack_secret_block: str = DEFAULT_SLACK_SECRET_BLOCK,
):
    """Run the assembly metadata update Nextflow pipeline in --screen_date mode: screens the
    metadata database for high-priority assemblies released after the specified date (YYYY-MM-DD),

    Output goes to {outdir}/{folder_name}; folder_name defaults to asm_update_screen_<screen_date>_<today>.
    """
    # Resolve the screen date and set the default output folder name.
    screen_date = datetime.strptime(screen_date, "%Y-%m-%d").strftime("%Y-%m-%d")
    default_folder = f"asm_update_screen_{screen_date}_{datetime.now().strftime('%Y-%m-%d')}"

    # Resolve credentials
    credentials = resolve_credentials(
        credentials, metadata_secret_block, slack_secret_block, slack_report, dry_run=dry_run
    )

    # Run the assembly metadata update pipeline for the specified screen date.
    run = asm_update_run(
        outdir=run_folder(outdir, folder_name, default_folder),
        enscode=resolve_enscode(enscode, dry_run),
        metadata_params_string=credentials.metadata_params_string,
        screen_date=screen_date,
        slack_report=slack_report,
        slack_params=credentials.slack_params,
    )
    return run_nextflow_pipeline(run, TaskRunOptions(dry_run=dry_run))
