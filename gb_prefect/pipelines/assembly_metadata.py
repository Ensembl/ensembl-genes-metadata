"""Spec for assembly_metadata pipeline."""

from datetime import datetime
from typing import Optional

from gb_prefect.models.nextflow import NextflowRun


def asm_registry_run(
    outdir: str,
    enscode: str,
    metadata_params_string: str,
    date: Optional[str] = None,
    gca_list: Optional[str] = None,
) -> NextflowRun:
    """Return the run of the assembly registry pipeline.

    Input modes:
    - gca_list: register exactly the accessions in that file, as given (--add_gca); no NCBI
      screening.
    - date (YYYY-MM-DD): screen NCBI for assemblies released after that date.
    - neither: screen NCBI from the DB's last regular update date.
    """
    return NextflowRun(
        name="assembly_metadata",
        enscode=enscode,
        pipeline="ensembl-genes-metadata/pipelines/assembly_metadata/main.nf",
        config="ensembl-genes-metadata/pipelines/assembly_metadata/nextflow.config",
        outdir=outdir,
        params={
            "output_dir": outdir,
            "enscode": enscode,
            "add_gca": True if gca_list else None,
            "gca_list": gca_list,
            "date": datetime.strptime(date, "%Y-%m-%d").strftime("%m/%d/%Y") if date else None,
        },
        secret_params={"metadata_params_string": metadata_params_string},
        with_report=True,
        dag_file="assembly_metadata_dag_%Y%m%d-%H%M%S.mmd",
    )
