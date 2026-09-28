import argparse
from pathlib import Path
from typing import Optional

from prefect import flow  # type: ignore

from gb_prefect.tasks.busco_gca import run_nextflow_busco_gca


@flow(name="BUSCO_gca", log_prints=True)
def busco_gca_flow(
    gca: str,
    taxon_id: str,
    outdir: str,
    busco_dataset: Optional[str] = None,
    enscode: Optional[str] = None,
    dry_run: bool = False,
    create_artifact: bool = True,
):
    """Run the BUSCO genome-statistics Nextflow pipeline for one GCA.

    This is the atomic unit deployed to the slurm-cli busco work pool -- one flow
    run submits as one Slurm job (see gb_prefect/worker/). Results are written
    under outdir/<gca> so concurrent runs against a shared base outdir don't
    collide.
    """
    run_outdir = str(Path(outdir) / gca)

    return run_nextflow_busco_gca(
        gca=gca,
        taxon_id=taxon_id,
        outdir=run_outdir,
        busco_dataset=busco_dataset,
        enscode=enscode,
        dry_run=dry_run,
        create_artifact=create_artifact,
    )


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
    args = parser.parse_args()

    busco_gca_flow(
        gca=args.gca,
        taxon_id=args.taxon_id,
        outdir=args.outdir,
        busco_dataset=args.busco_dataset,
        enscode=args.enscode,
        dry_run=args.dry_run,
        create_artifact=args.create_artifact,
    )
