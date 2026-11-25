"""Write and validate planner-generated RNA-seq evidence job manifests."""

import json
import os
import tempfile
from pathlib import Path
from typing import Iterable, List

from gb_prefect.models.evidence import EvidenceJob


def write_job_manifest(path: str, jobs: Iterable[EvidenceJob]) -> str:
    """Atomically write one validated JSON object per evidence job.

    The planner owns discovery and must pass explicit paths. This helper only validates the job
    contract and writes a deterministic manifest; it never scans a production directory.
    """
    target = Path(path)
    if not target.is_absolute():
        raise ValueError(f"job manifest must be an absolute path: {path}")
    records: List[dict] = []
    seen = set()
    for job in jobs:
        if job.job_id in seen:
            raise ValueError(f"duplicate evidence job_id: {job.job_id}")
        seen.add(job.job_id)
        records.append(job.to_mapping())
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=".evidence_jobs.", dir=str(target.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            for record in records:
                handle.write(json.dumps(record, sort_keys=True) + "\n")
        os.chmod(temporary, 0o600)
        os.replace(temporary, target)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    return str(target)
