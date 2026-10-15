"""Spec for ensembl-genes-nf's statistics pipeline in BUSCO genome mode (--run_busco_ncbi)."""

from typing import Any, Mapping, Optional

from gb_prefect.models.nextflow import CsvInput, NextflowRun


def busco_run(  # pylint: disable=too-many-arguments
    gca: str,
    taxon_id: str,
    outdir: str,
    enscode: str,
    busco_dataset: Optional[str] = None,
    update_registry: bool = False,
    asm_metadata: Optional[Mapping[str, Any]] = None,
) -> NextflowRun:
    """Return the run of the statistics pipeline's NCBI BUSCO mode (--run_busco_ncbi) for one GCA.

    outdir is this GCA's own run directory. busco_dataset picks the BUSCO lineage; omitted, the
    pipeline chooses one from the taxon.

    update_registry=True adds --update_registry: the pipeline loads the results into the
    assembly metadata DB and marks genome_busco.status as done. It needs asm_metadata (DB
    connection params), which goes to the owner-only -params-file -- never the command line.
    """
    if update_registry and asm_metadata is None:
        raise ValueError("asm_metadata is required when update_registry is True")

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
        secret_params={"asm_metadata": asm_metadata} if update_registry else {},
    )
