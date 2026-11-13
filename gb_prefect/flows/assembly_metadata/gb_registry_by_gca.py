from datetime import datetime
from typing import Optional

from prefect import flow  # type: ignore

from gb_prefect.models.pipeline_options import TaskRunOptions
from gb_prefect.pipelines.assembly_metadata import asm_registry_run
from gb_prefect.tasks.nextflow import run_nextflow_pipeline
from gb_prefect.utils.credentials_utils import PipelineCredentials, resolve_credentials
from gb_prefect.utils.enscode_utils import resolve_enscode
from gb_prefect.utils.path_utils import check_input_file, run_folder


@flow(name="gb_registry_by_gca", log_prints=True)
def gb_registry_by_gca_flow(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    gca_list: str,
    outdir: str,
    folder_name: Optional[str] = None,
    enscode: Optional[str] = None,
    dry_run: bool = False,
    credentials: Optional[PipelineCredentials] = None,
    metadata_secret_block: Optional[str] = None,
):
    """
    Run the assembly registry Nextflow pipeline for the GCA accessions listed in gca_list
    (--add_gca mode, the whole file as given).

    gca_list must be an absolute path. Output goes to {outdir}/{folder_name}; folder_name
    defaults to asm_registry_gca_<today>.
    """
    # Resolve credentials
    credentials = resolve_credentials(
        credentials, dry_run=dry_run, metadata_secret_block=metadata_secret_block
    )

    # Run the assembly registry pipeline for the specified GCA accessions.
    run = asm_registry_run(
        outdir=run_folder(outdir, folder_name, f"asm_registry_gca_{datetime.now().strftime('%Y-%m-%d')}"),
        enscode=resolve_enscode(enscode, dry_run),
        metadata_params_string=credentials.metadata_params_string,
        gca_list=check_input_file(gca_list, dry_run),
    )
    return run_nextflow_pipeline(run, TaskRunOptions(dry_run=dry_run))
