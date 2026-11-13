from pathlib import Path
from typing import Optional

from prefect import flow  # type: ignore

from gb_prefect.models.pipeline_options import TaskRunOptions
from gb_prefect.pipelines.statistics_busco_genome import busco_run
from gb_prefect.tasks.statistics_cleanup import cleanup_statistics_outdir, record_busco_failure
from gb_prefect.tasks.nextflow import run_nextflow_pipeline
from gb_prefect.utils.credentials_utils import PipelineCredentials, resolve_credentials
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
    metadata_secret_block: Optional[str] = None,
    post_run_cleanup: bool = False,
):
    """Run the statistics Nextflow pipeline's BUSCO genome mode for one GCA, under <outdir>/<gca>.

    This is the unit the BUSCO dispatchers trigger, one flow run (one Slurm job) per GCA; the
    per-GCA folder keeps concurrent runs against the same outdir apart.

    gca, taxon_id: the assembly and its NCBI taxon. busco_dataset: the BUSCO lineage; omitted,
        the pipeline chooses one from the taxon.
    enscode: ENSCODE directory holding ensembl-genes-nf; defaults to $ENSCODE.
    dry_run: write and log the command without running Nextflow or loading credentials.
    create_artifact: publish the run-summary artifact.
    update_registry: load the results into the assembly metadata DB and mark the GCA's
        genome_busco status done (the pipeline's --update_registry). Needs the DB credentials:
        credentials, or the metadata_secret_block Secret block (the standard one if None).
    post_run_cleanup: set by the automatic dispatcher. On success, remove <outdir>/<gca>; on any
        failure (including failing to load credentials), mark the GCA's status failed, remove
        Nextflow's scratch and the credentials file but keep the logs, and re-raise so the run
        shows Failed. A Slurm job killed outright (time limit, out of memory) never reaches this
        step, so its GCA stays in_progress.
    """
    run_outdir = str(Path(outdir) / gca)
    cleanup = post_run_cleanup and not dry_run

    try:
        # DB credentials only for --update_registry; inside the try, so a failure to load them
        # still marks the GCA failed. On a dry run they are a placeholder (nothing loaded).
        metadata_params_string = None
        if update_registry:
            metadata_params_string = resolve_credentials(
                credentials, dry_run=dry_run, metadata_secret_block=metadata_secret_block
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
