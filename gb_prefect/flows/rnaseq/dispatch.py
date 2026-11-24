"""Capacity-aware dispatcher for planner-generated RNA-seq evidence jobs."""

import json
import time
from pathlib import Path
from typing import Dict, List, Optional

import httpx
from prefect import flow  # type: ignore
from prefect.deployments import run_deployment  # type: ignore
from prefect.exceptions import ObjectNotFound  # type: ignore
from prefect.runtime import flow_run  # type: ignore

from gb_prefect.models.evidence import (
    EvidenceKind,
    EvidenceJob,
    EvidenceStatus,
    read_job_manifest,
    transition_state,
)
from gb_prefect.tasks.busco_dispatch import get_free_slots


DEFAULT_TRACK_DEPLOYMENT = "RNAseq_track_preparation/RNAseq_track_preparation"
DEFAULT_ENA_DEPLOYMENT = "RNAseq_ena_submission/RNAseq_ena_submission"
TRIGGER_ATTEMPTS = 3
TRIGGER_RETRY_WAIT_SECONDS = 5


class TriggerError(RuntimeError):
    """A child deployment trigger failed, with whether the child may already exist."""

    def __init__(self, error: Exception, run_may_exist: bool):
        super().__init__(str(error))
        self.run_may_exist = run_may_exist


def _job_file(manifest: str, job: EvidenceJob) -> str:
    path = Path(manifest).parent / "jobs" / f"{job.job_id}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        existing = EvidenceJob.from_mapping(json.loads(path.read_text(encoding="utf-8")))
        if existing.to_mapping() != job.to_mapping():
            raise ValueError(f"Existing job file does not match manifest for {job.job_id}: {path}")
    else:
        path.write_text(json.dumps(job.to_mapping(), indent=2, sort_keys=True) + "\n", encoding="utf-8")
        path.chmod(0o600)
    return str(path)


def _trigger(
    deployment: str,
    job_file: str,
    job: EvidenceJob,
    idempotency_key: str,
    webin_user: Optional[str],
) -> str:
    """Trigger a child deployment without blocking the dispatcher Slurm/process slot."""
    parameters = {"job_file": job_file}
    if job.kind is EvidenceKind.ENA:
        if not webin_user:
            raise ValueError("webin_user is required when dispatching ENA jobs")
        parameters["webin_user"] = webin_user
    error: Optional[Exception] = None
    for attempt in range(1, TRIGGER_ATTEMPTS + 1):
        try:
            child = run_deployment(
                name=deployment,
                parameters=parameters,
                timeout=0,
                idempotency_key=idempotency_key,
            )
            return str(getattr(child, "id", None) or "unknown")  # type: ignore[union-attr]
        except ObjectNotFound as err:
            raise TriggerError(err, run_may_exist=False) from err
        except httpx.HTTPStatusError as err:
            if 400 <= err.response.status_code < 500 and err.response.status_code not in (408, 429):
                raise TriggerError(err, run_may_exist=False) from err
            error = err
        except Exception as err:  # pylint: disable=broad-exception-caught
            error = err
        if attempt < TRIGGER_ATTEMPTS:
            print(f"Trigger attempt {attempt}/{TRIGGER_ATTEMPTS} failed: {error!r}")
            time.sleep(TRIGGER_RETRY_WAIT_SECONDS)
    assert error is not None
    raise TriggerError(error, run_may_exist=True) from error


@flow(name="RNAseq_evidence_dispatch", log_prints=True)
# pylint: disable=too-many-locals
def rnaseq_evidence_dispatch_flow(
    job_manifest: str,
    track_deployment: str = DEFAULT_TRACK_DEPLOYMENT,
    ena_deployment: str = DEFAULT_ENA_DEPLOYMENT,
    dry_run: bool = False,
    max_jobs_per_kind: Optional[int] = None,
    webin_user: Optional[str] = None,
) -> Dict[str, List[str]]:
    """Dispatch planned JSONL evidence jobs up to each child deployment's free capacity.

    The dispatcher never discovers files and never submits Slurm jobs itself. It claims a job's
    shared state file before triggering the matching child deployment; the child flow owns the
    Nextflow run and moves the state to running/succeeded/failed.
    """
    jobs = read_job_manifest(job_manifest)
    if not dry_run and any(job.kind is EvidenceKind.ENA for job in jobs) and not webin_user:
        raise ValueError("webin_user is required when dispatching ENA jobs")
    deployments = {EvidenceKind.TRACK: track_deployment, EvidenceKind.ENA: ena_deployment}
    result: Dict[str, List[str]] = {"triggered": [], "skipped": [], "failed": []}
    dispatcher_flow_id = getattr(flow_run, "id", None)
    dispatcher_id = str(dispatcher_flow_id) if dispatcher_flow_id else "direct"

    for kind, deployment in deployments.items():
        candidates = [job for job in jobs if job.kind is kind]
        if not candidates:
            continue
        slots = max_jobs_per_kind if max_jobs_per_kind is not None else get_free_slots(deployment)
        if slots <= 0:
            result["skipped"].extend(job.job_id for job in candidates)
            continue
        triggered_count = 0
        for job in candidates:
            if triggered_count >= slots:
                result["skipped"].append(job.job_id)
                continue
            state = EvidenceStatus(read_job_manifest_state(job))
            if state not in (EvidenceStatus.PLANNED, EvidenceStatus.FAILED):
                result["skipped"].append(job.job_id)
                continue
            if dry_run:
                result["triggered"].append(f"{kind.value}:{job.job_id}")
                continue
            job_file = _job_file(job_manifest, job)
            transition_state(job, EvidenceStatus.DISPATCHED, flow_run_id=dispatcher_id)
            try:
                child_id = _trigger(deployment, job_file, job, f"{dispatcher_id}-{job.job_id}", webin_user)
            except TriggerError as err:
                if not err.run_may_exist:
                    transition_state(job, EvidenceStatus.FAILED, flow_run_id=dispatcher_id, error=str(err))
                result["failed"].append(job.job_id)
                continue
            result["triggered"].append(f"{kind.value}:{job.job_id}:{child_id}")
            triggered_count += 1
    return result


def read_job_manifest_state(job: EvidenceJob) -> str:
    """Read a state value without importing Prefect task machinery into the dispatcher loop."""
    if not job.state_path.exists():
        return EvidenceStatus.PLANNED.value
    value = json.loads(job.state_path.read_text(encoding="utf-8"))
    return value.get("status", EvidenceStatus.PLANNED.value)
