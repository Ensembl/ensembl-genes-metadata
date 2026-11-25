"""Nextflow specs for RNA-seq browser-track preparation and ENA submission."""

from typing import Optional

from gb_prefect.models.evidence import EvidenceJob
from gb_prefect.models.nextflow import NextflowRun


def track_preparation_run(job: EvidenceJob, enscode: str) -> NextflowRun:
    """Build one track-preparation run for an explicit GCA/entity manifest."""
    outdir = str(job.outdir / "track")
    params = {
        "input_manifest": job.track_manifest,
        "chrom_sizes": job.chrom_sizes,
        "assembly_release": job.assembly_release,
        "outdir": outdir,
        "track_types": job.track_types,
        "manifest_base_dir": job.manifest_base_dir,
        "lts_stage_dir": job.lts_stage_dir,
    }
    params.update(job.extra_params)
    return NextflowRun(
        name=f"rnaseq_track_{job.gca_accession}",
        enscode=enscode,
        pipeline="ensembl-genes-nf/pipelines/track_file_preparation/main.nf",
        config="ensembl-genes-nf/pipelines/track_file_preparation/nextflow.config",
        outdir=outdir,
        params={key: value for key, value in params.items() if value is not None},
        profile=job.nextflow_profile,
        with_report=True,
        dag_file="rnaseq_track_dag_%Y%m%d-%H%M%S.mmd",
    )


def ena_submission_run(
    job: EvidenceJob,
    enscode: str,
    webin_user: str,
    webin_password: Optional[str] = None,
) -> NextflowRun:
    """Build one ENA submission run for an explicit project/file manifest.

    ``webin_password`` is written to the private Nextflow secrets store by the generic runner;
    it never becomes a command-line or params-file value. The ENA pipeline's ``secret`` directive
    reads the resulting ``ENA_WEBIN_PASSWORD`` key.
    """
    outdir = str(job.outdir / "ena")
    params = {
        "manifest": job.ena_manifest,
        "mode": job.mode.value,
        "webin_user": webin_user,
        "outdir": outdir,
        "remote_dir": job.remote_dir,
        "upload_parallelism": job.upload_parallelism,
        "reheader_bams": job.reheader_bams,
        "convert_to_cram": job.convert_to_cram,
        "create_project": job.create_project,
        "reference_fasta": job.reference_fasta,
        "reference_assembly_report": job.reference_assembly_report,
    }
    params.update(job.extra_params)
    return NextflowRun(
        name=f"rnaseq_ena_{job.gca_accession}",
        enscode=enscode,
        pipeline="ensembl-genes-nf/pipelines/ena_submit/main.nf",
        config="ensembl-genes-nf/pipelines/ena_submit/nextflow.config",
        outdir=outdir,
        params={key: value for key, value in params.items() if value is not None},
        profile=job.nextflow_profile,
        secret_store={"ENA_WEBIN_PASSWORD": webin_password} if webin_password else {},
        with_report=True,
        dag_file="rnaseq_ena_dag_%Y%m%d-%H%M%S.mmd",
    )
