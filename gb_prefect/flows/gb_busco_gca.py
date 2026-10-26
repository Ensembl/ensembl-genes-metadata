from pathlib import Path
from typing import Optional

from prefect import flow  # type: ignore

from gb_prefect.models.pipeline_options import PipelineCredentials, TaskRunOptions
from gb_prefect.pipelines.statistics_busco_genome import busco_run
from gb_prefect.tasks.statistics_cleanup import cleanup_statistics_outdir, record_busco_failure
from gb_prefect.tasks.nextflow import run_nextflow_pipeline
from gb_prefect.utils.credentials_utils import (
    DEFAULT_METADATA_SECRET_BLOCK,
    DEFAULT_SLACK_SECRET_BLOCK,
    resolve_credentials,
)
from gb_prefect.utils.enscode_utils import resolve_enscode


@flow(name="BUSCO_gca", log_prints=True)
def busco_gca_flow(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    gca: str,
    taxon_id: str,
    outdir: str,
    busco_dataset: Optional[str] = None,
    enscode: Optional[str] = None,
    dry_run: bool = False,
    create_artifact: bool = True,
    update_registry: bool = False,
    credentials: Optional[PipelineCredentials] = None,
    metadata_secret_block: str = DEFAULT_METADATA_SECRET_BLOCK,
    post_run_cleanup: bool = False,
):
    """
    Run the BUSCO genome from statistics Nextflow pipeline for one GCA.
    """
    run_outdir = str(Path(outdir) / gca)
    cleanup = post_run_cleanup and not dry_run

    try:
        # DB credentials only for --update_registry; inside the try, so a failure to load them
        # still marks the GCA failed. On a dry run they are a placeholder (nothing loaded).
        metadata_params_string = None
        if update_registry:
            metadata_params_string = resolve_credentials(
                credentials,
                metadata_secret_block,
                DEFAULT_SLACK_SECRET_BLOCK,
                slack_report=False,
                dry_run=dry_run,
            ).metadata_params_string

        run = busco_run(
            gca=gca,
            taxon_id=taxon_id,
            outdir=run_outdir,
            enscode=resolve_enscode(enscode, dry_run),
            busco_dataset=busco_dataset,
            update_registry=update_registry,
            metadata_params_string=metadata_params_string,
        )
        result = run_nextflow_pipeline(run, TaskRunOptions(dry_run=dry_run, create_artifact=create_artifact))
    except Exception:
        if cleanup:
            try:
                record_busco_failure(gca, credentials, metadata_secret_block)
            finally:
                cleanup_statistics_outdir(run_outdir, success=False)
        raise

    if cleanup:
        cleanup_statistics_outdir(run_outdir, success=True)
    return result
