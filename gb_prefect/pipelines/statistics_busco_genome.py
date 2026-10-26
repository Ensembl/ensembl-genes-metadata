"""Spec for ensembl-genes-nf's statistics pipeline in BUSCO genome mode (--run_busco_ncbi)."""

from typing import Optional

from gb_prefect.models.nextflow import CsvInput, NextflowRun


def busco_run(  # pylint: disable=too-many-arguments
    gca: str,
    taxon_id: str,
    outdir: str,
    enscode: str,
    busco_dataset: Optional[str] = None,
    update_registry: bool = False,
    metadata_params_string: Optional[str] = None,
) -> NextflowRun:
    """Return the run of the statistics pipeline's NCBI BUSCO mode (--run_busco_ncbi) for one GCA.

    outdir is this GCA's own run directory. busco_dataset picks the BUSCO lineage; omitted, the
    pipeline chooses one from the taxon.

    update_registry=True adds --update_registry: the pipeline loads the results into the
    assembly metadata DB and marks genome_busco.status as done. It needs metadata_params_string
    (DB connection params, JSON), passed as the pipeline's asm_metadata in the owner-only
    -params-file -- never the command line. The pipeline accepts asm_metadata as a JSON string
    or object.
    """
    if update_registry and metadata_params_string is None:
        raise ValueError("metadata_params_string is required when update_registry is True")

    return NextflowRun(
        name=gca,
        enscode=enscode,
        pipeline="ensembl-genes-nf/pipelines/statistics/main.nf",
        outdir=outdir,
        csv_input=CsvInput(
            param="csvFile",
            columns=("gca", "taxon_id", "busco_dataset"),
            rows=[{"gca": gca, "taxon_id": taxon_id, "busco_dataset": busco_dataset}],
        ),
        params={
            "run_busco_ncbi": True,
            "outdir": outdir,
            "enscode": enscode,
            "buscoDataset": busco_dataset or None,
            "update_registry": update_registry or None,
        },
        secret_params={"asm_metadata": metadata_params_string} if update_registry else {},
        metro_map="ensembl-genes-nf/pipelines/statistics/assets/metro_map_busco_genome.mmd",
    )
