"""Spec for assembly metadata update pipeline ."""

from datetime import datetime
from typing import Optional

from gb_prefect.models.nextflow import NextflowRun


def asm_update_run(  # pylint: disable=too-many-arguments
    outdir: str,
    enscode: str,
    metadata_params_string: str,
    gca_list: Optional[str] = None,
    screen_date: Optional[str] = None,
    slack_report: bool = True,
    slack_params: Optional[str] = None,
) -> NextflowRun:
    """Return the NextflowRun of the assembly metadata update pipeline.

    Exactly one input mode:
    - gca_list: check/update the accessions in that file, as given (--gca_input).
    - screen_date (YYYY-MM-DD): screen the metadata DB for high-priority assemblies released
      after that date.
    """
    # Validate input arguments and prepare secret parameters.
    if bool(gca_list) == bool(screen_date):
        raise ValueError("Give exactly one of gca_list or screen_date")
    if slack_report and not slack_params:
        raise ValueError("slack_params is required when slack_report is True")
    if screen_date:
        screen_date = datetime.strptime(screen_date, "%Y-%m-%d").strftime("%Y-%m-%d")
    secret_params = {"metadata_params_string": metadata_params_string}
    if slack_report and slack_params:
        secret_params["slack_params"] = slack_params

    # Return the configured NextflowRun object.
    return NextflowRun(
        name="assembly_metadata_update",
        enscode=enscode,
        pipeline="ensembl-genes-metadata/pipelines/assembly_metadata_update/main.nf",
        config="ensembl-genes-metadata/pipelines/assembly_metadata_update/nextflow.config",
        outdir=outdir,
        params={
            "output_dir": outdir,
            "gca_input": True if gca_list else None,
            "gca_list": gca_list,
            "screen_date": screen_date,
            "slack_report": slack_report,
        },
        secret_params=secret_params,
        with_report=True,
        dag_file="assembly_metadata_update_dag_%Y%m%d-%H%M%S.mmd",
        metro_map="ensembl-genes-metadata/pipelines/assembly_metadata_update/assets/metro_map.mmd",
    )
