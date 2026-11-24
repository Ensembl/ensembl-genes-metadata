"""Prefect deployment entry point for one RNA-seq track job."""

import json
from pathlib import Path
from typing import Optional

from prefect import flow  # type: ignore
from prefect.runtime import flow_run  # type: ignore

from gb_prefect.models.evidence import EvidenceKind, EvidenceJob, EvidenceStatus, transition_state
from gb_prefect.models.pipeline_options import TaskRunOptions
from gb_prefect.pipelines.rnaseq_evidence import track_preparation_run
from gb_prefect.tasks.nextflow import run_nextflow_pipeline
from gb_prefect.utils.enscode_utils import resolve_enscode


def _load_job(job_file: str) -> EvidenceJob:
    value = json.loads(Path(job_file).read_text(encoding="utf-8"))
    job = EvidenceJob.from_mapping(value)
    if job.kind is not EvidenceKind.TRACK:
        raise ValueError(f"Job {job.job_id} is {job.kind.value}, expected track")
    return job


@flow(name="RNAseq_track_preparation", log_prints=True, retries=1, retry_delay_seconds=60)
def rnaseq_track_preparation_flow(
    job_file: str,
    enscode: Optional[str] = None,
    dry_run: bool = False,
    create_artifact: bool = True,
):
    """Run track preparation for one planner-generated JSON job record.

    The job record and output directory are stable across Prefect retries, so the generic runner
    adds Nextflow ``-resume`` and reuses completed entity work.
    """
    job = _load_job(job_file)
    flow_id = str(flow_run.id) if flow_run.id else None
    if not dry_run:
        transition_state(job, EvidenceStatus.RUNNING, flow_run_id=flow_id)
    try:
        run = track_preparation_run(job, resolve_enscode(enscode, dry_run))
        result = run_nextflow_pipeline(run, TaskRunOptions(dry_run=dry_run, create_artifact=create_artifact))
    except Exception as err:
        if not dry_run:
            transition_state(job, EvidenceStatus.FAILED, flow_run_id=flow_id, error=str(err))
        raise
    if not dry_run:
        transition_state(job, EvidenceStatus.SUCCEEDED, flow_run_id=flow_id)
    return result
