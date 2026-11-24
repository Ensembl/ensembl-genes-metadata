import json
import stat
from types import SimpleNamespace
from pathlib import Path

import pytest

from gb_prefect.flows.rnaseq.dispatch import rnaseq_evidence_dispatch_flow
from gb_prefect.models.evidence import (
    EvidenceJob,
    EvidenceKind,
    EvidenceStatus,
    read_job_manifest,
    read_state,
    transition_state,
)
from gb_prefect.nextflow.command import nextflow_env, write_secret_store
from gb_prefect.pipelines.rnaseq_evidence import ena_submission_run, track_preparation_run


def _job(tmp_path: Path, kind: EvidenceKind, job_id: str) -> EvidenceJob:
    common = {
        "job_id": job_id,
        "kind": kind,
        "gca_accession": "GCA_000001405.28",
        "assembly_release": "GCA_000001405.28-Ensembl-test",
        "output_root": str(tmp_path / "outputs"),
    }
    if kind is EvidenceKind.TRACK:
        common.update(
            track_manifest=str(tmp_path / "track.tsv"),
            chrom_sizes=str(tmp_path / "chrom.sizes"),
        )
    else:
        common.update(ena_manifest=str(tmp_path / "ena.tsv"))
    return EvidenceJob(**common)


def test_job_manifest_and_state_transitions_are_idempotent(tmp_path):
    job = _job(tmp_path, EvidenceKind.TRACK, "track-GCA_000001405.28")
    manifest = tmp_path / "jobs.jsonl"
    manifest.write_text(json.dumps(job.to_mapping()) + "\n", encoding="utf-8")

    [loaded] = read_job_manifest(str(manifest))
    assert loaded == job
    assert read_state(job)["status"] == "planned"
    transition_state(job, EvidenceStatus.DISPATCHED, flow_run_id="dispatcher")
    transition_state(job, EvidenceStatus.RUNNING, flow_run_id="child")
    final = transition_state(job, EvidenceStatus.SUCCEEDED, flow_run_id="child")
    assert final["status"] == "succeeded"
    assert [item["status"] for item in final["history"]] == ["dispatched", "running", "succeeded"]

    with pytest.raises(ValueError, match="Invalid evidence state transition"):
        transition_state(job, EvidenceStatus.RUNNING)


def test_specs_keep_pipeline_fanout_inside_nextflow(tmp_path):
    track = _job(tmp_path, EvidenceKind.TRACK, "track-job")
    track_run = track_preparation_run(track, "/hps/enscode")
    assert track_run.pipeline.endswith("track_file_preparation/main.nf")
    assert track_run.params["input_manifest"] == track.track_manifest
    assert track_run.outdir.endswith("/track-job/track")

    ena = _job(tmp_path, EvidenceKind.ENA, "ena-job")
    ena_run = ena_submission_run(ena, "/hps/enscode", "Webin-123", "p@ssword")
    assert ena_run.pipeline.endswith("ena_submit/main.nf")
    assert ena_run.params["manifest"] == ena.ena_manifest
    assert ena_run.profile is None
    assert ena_run.secret_store == {"ENA_WEBIN_PASSWORD": "p@ssword"}
    assert nextflow_env(ena_run)["NXF_HOME"].endswith("/ena-job/ena/.nextflow_home")


def test_secret_store_is_private_and_json_compatible(tmp_path):
    path = tmp_path / ".nextflow_home" / "secrets" / "store.json"
    write_secret_store(path, {"ENA_WEBIN_PASSWORD": "p@ssword"})
    assert stat.S_IMODE(path.stat().st_mode) == 0o600
    assert json.loads(path.read_text(encoding="utf-8")) == [
        {"name": "ENA_WEBIN_PASSWORD", "value": "p@ssword"}
    ]


def test_dispatch_dry_run_respects_per_kind_capacity(tmp_path):
    track = _job(tmp_path, EvidenceKind.TRACK, "track-job")
    ena = _job(tmp_path, EvidenceKind.ENA, "ena-job")
    manifest = tmp_path / "jobs.jsonl"
    manifest.write_text(
        "\n".join(json.dumps(job.to_mapping()) for job in (track, ena)) + "\n", encoding="utf-8"
    )

    result = rnaseq_evidence_dispatch_flow.fn(str(manifest), dry_run=True, max_jobs_per_kind=1)
    assert sorted(item.split(":", 2)[:2] for item in result["triggered"]) == [
        ["ena", "ena-job"],
        ["track", "track-job"],
    ]


def test_dispatch_forwards_webin_user_and_claims_job(tmp_path, monkeypatch):
    ena = _job(tmp_path, EvidenceKind.ENA, "ena-job")
    manifest = tmp_path / "jobs.jsonl"
    manifest.write_text(json.dumps(ena.to_mapping()) + "\n", encoding="utf-8")
    calls = []

    def fake_run_deployment(**kwargs):
        calls.append(kwargs)
        return SimpleNamespace(id="child-run")

    monkeypatch.setattr("gb_prefect.flows.rnaseq.dispatch.run_deployment", fake_run_deployment)
    result = rnaseq_evidence_dispatch_flow.fn(
        str(manifest), dry_run=False, max_jobs_per_kind=1, webin_user="Webin-123"
    )

    assert result["triggered"] == ["ena:ena-job:child-run"]
    assert calls[0]["parameters"]["webin_user"] == "Webin-123"
    assert read_state(ena)["status"] == "dispatched"


def test_real_dispatch_requires_webin_user_for_ena_jobs(tmp_path):
    ena = _job(tmp_path, EvidenceKind.ENA, "ena-job")
    manifest = tmp_path / "jobs.jsonl"
    manifest.write_text(json.dumps(ena.to_mapping()) + "\n", encoding="utf-8")

    with pytest.raises(ValueError, match="webin_user is required"):
        rnaseq_evidence_dispatch_flow.fn(str(manifest), dry_run=False, max_jobs_per_kind=1)
