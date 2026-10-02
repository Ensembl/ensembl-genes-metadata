import argparse
from pathlib import Path
from typing import Optional

from prefect import flow  # type: ignore

from gb_prefect.models.pipeline_options import PipelineCredentials
from gb_prefect.tasks.busco_gca import (
    cleanup_busco_outdir,
    record_busco_failure,
    run_nextflow_busco_gca,
)
from gb_prefect.utils.credentials_utils import DEFAULT_METADATA_SECRET_BLOCK


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
    """Run the BUSCO genome-statistics Nextflow pipeline for one GCA.

    This is the atomic unit deployed to the slurm-cli busco work pool -- one flow
    run submits as one Slurm job (see gb_prefect/worker/). Results are written
    under outdir/<gca> so concurrent runs against a shared base outdir don't
    collide.

    update_registry=True passes --update_registry to the pipeline, which loads the genome
    results into the assembly metadata DB and marks genome_busco.status as done on
    completion. DB connection params come from credentials, or from the
    metadata_secret_block Prefect Secret block when omitted (the deployment case) -- see
    gb_prefect.tasks.busco_gca.run_nextflow_busco_gca.

    post_run_cleanup=True (always set by the automatic dispatcher, optional in the others) adds a step
    after the pipeline: on success the whole <outdir>/<gca> directory is removed; on failure
    genome_busco.status is set to failed, Nextflow scratch and the credentials file are
    removed (logs kept), and the original error is re-raised so the run still shows Failed.
    A Slurm job killed outright (time limit, OOM) never reaches this step, so its GCA stays
    in_progress.
    """
    run_outdir = str(Path(outdir) / gca)

    run_kwargs = {
        "gca": gca,
        "taxon_id": taxon_id,
        "outdir": run_outdir,
        "busco_dataset": busco_dataset,
        "enscode": enscode,
        "dry_run": dry_run,
        "create_artifact": create_artifact,
        "update_registry": update_registry,
        "credentials": credentials,
        "metadata_secret_block": metadata_secret_block,
    }
    if not post_run_cleanup or dry_run:
        return run_nextflow_busco_gca(**run_kwargs)

    try:
        result = run_nextflow_busco_gca(**run_kwargs)
    except Exception:
        try:
            record_busco_failure(gca, credentials, metadata_secret_block)
        finally:
            cleanup_busco_outdir(run_outdir, success=False)
        raise

    cleanup_busco_outdir(run_outdir, success=True)
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
