import argparse
from pathlib import Path
from typing import Optional

from prefect import flow  # type: ignore

from gb_prefect.models.pipeline_options import PipelineCredentials, TaskRunOptions
from gb_prefect.pipelines.statistics_busco_genome import busco_run
from gb_prefect.tasks.statistics_cleanup import cleanup_statistics_outdir, record_busco_failure
from gb_prefect.tasks.nextflow import run_nextflow_pipeline
from gb_prefect.utils.credentials_utils import DEFAULT_METADATA_SECRET_BLOCK, load_metadata_db_params
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
        asm_metadata = {}
        if update_registry and not dry_run:
            asm_metadata = load_metadata_db_params(credentials, metadata_secret_block)

        run = busco_run(
            gca=gca,
            taxon_id=taxon_id,
            outdir=run_outdir,
            enscode=resolve_enscode(enscode, dry_run),
            busco_dataset=busco_dataset,
            update_registry=update_registry,
            asm_metadata=asm_metadata,
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--gca", required=True, help="GCA accession to run BUSCO for.")
    parser.add_argument("--taxon-id", required=True, help="NCBI taxon ID for the GCA.")
    parser.add_argument("--outdir", required=True, help="Base output directory.")
    parser.add_argument(
        "--busco-dataset", help="BUSCO lineage dataset to use; omit to let the pipeline choose."
    )
    parser.add_argument("--enscode", help="Path to the ENSCODE directory.")
    parser.add_argument(
        "--dry-run", action="store_true", help="Create the Nextflow command without running it."
    )
    parser.add_argument(
        "--create-artifact", action="store_true", help="Create artifact after running the pipeline."
    )
    parser.add_argument(
        "--update_registry",
        action="store_true",
        help="Load results into the assembly metadata DB and mark genome_busco.status as done.",
    )
    parser.add_argument(
        "--metadata-params-string",
        required=False,
        help="JSON string with metadata DB connection parameters, used with --update_registry. "
        "If omitted, loaded from the Prefect Secret block.",
    )
    parser.add_argument(
        "--post-run-cleanup",
        action="store_true",
        help="After the run: on success remove <outdir>/<gca>; on failure mark the GCA as "
        "failed and remove Nextflow scratch (logs kept). Set by the automatic dispatcher.",
    )
    args = parser.parse_args()

    busco_gca_flow(
        gca=args.gca,
        taxon_id=args.taxon_id,
        outdir=args.outdir,
        busco_dataset=args.busco_dataset,
        enscode=args.enscode,
        dry_run=args.dry_run,
        create_artifact=args.create_artifact,
        update_registry=args.update_registry,
        credentials=(
            PipelineCredentials(metadata_params_string=args.metadata_params_string)
            if args.metadata_params_string
            else None
        ),
        post_run_cleanup=args.post_run_cleanup,
    )
