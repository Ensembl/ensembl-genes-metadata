# RNA-seq evidence orchestration

This layer orchestrates the two pipelines in `ensembl-genes-nf` without moving entity/BAM
fan-out into Prefect:

```text
planner → JSONL EvidenceJob manifest on shared storage
        → Prefect dispatcher claims jobs and respects pool capacity
        → one Slurm flow run per GCA/release workload
        → Nextflow fans out over track entities or ENA files
        → evidence_state.json + pipeline outputs + reports
```

The schedulable unit is a GCA plus Ensembl release, not an individual BAM. Track and ENA jobs use
separate work pools because their bottlenecks differ. A planner may emit two records for the same
release, one with `kind: "track"` and one with `kind: "ena"`, using distinct `job_id` values.

## Job contract

The planner writes one JSON object per line. Paths are absolute and explicit; the dispatcher never
scans a production directory.

```json
{
  "job_id": "track_GCA_000001405.28_Ensembl_2026-10",
  "kind": "track",
  "gca_accession": "GCA_000001405.28",
  "assembly_release": "GCA_000001405.28-Ensembl-2026-10",
  "output_root": "/hps/.../rnaseq-evidence",
  "track_manifest": "/hps/.../track_manifest.tsv",
  "chrom_sizes": "/hps/.../chrom.sizes",
  "track_types": "all",
  "nextflow_profile": "singularity,slurm",
  "lts_stage_dir": "/hps/.../staged"
}
```

An ENA record uses `kind: "ena"` and `ena_manifest` instead of `track_manifest`/`chrom_sizes`,
plus `mode`, `remote_dir`, `upload_parallelism`, and the reference/reheader policy. Credentials
are never present in this file.

Each job gets a stable `<output_root>/<job_id>/evidence_state.json`. State transitions are locked
and atomic, so dispatch retries cannot claim the same job simultaneously. Prefect retries reuse
the same job directory and the generic runner adds Nextflow `-resume`.

## Deployments and pools

Deploy the following entries from `gb_prefect/setup/deploy_flows.py`:

- `RNAseq_track_preparation` on `genebuild-rnaseq-track` (`slurm-cli`);
- `RNAseq_ena_submission` on `genebuild-rnaseq-ena` (`slurm-cli`);
- `RNAseq_evidence_dispatch` on `genebuild-dispatch` (`process`).

The Slurm flow job is the Nextflow head job. The pipeline's process resources and executor remain
in `ensembl-genes-nf`; the head-job time limit must cover the whole Nextflow run. Set the ENA pool
limit and per-run `upload_parallelism` together so the product of concurrent runs and uploads stays
within the Webin/FTP operating limit.

## Secrets

Create the Webin Secret block interactively:

```bash
python gb_prefect/setup/create_secrets.py --webin-password
```

The ENA flow writes `ENA_WEBIN_PASSWORD` into a private `NXF_HOME/secrets/store.json` with mode
0600. The password is not a deployment parameter, command-line argument, command record, log, or
Prefect artifact. The source ENA pipeline must continue to use Nextflow's `secret
'ENA_WEBIN_PASSWORD'` directive.

## Dry-run and dispatch

Create a manifest with `gb_prefect.models.evidence.EvidenceJob` and
`gb_prefect.utils.evidence_manifest.write_job_manifest`, then preview a child flow:

```bash
python -c "from gb_prefect.flows.rnaseq.ena_submission import rnaseq_ena_submission_flow as f; f(job_file='/abs/jobs/ena.json', webin_user='Webin-XXXXX', dry_run=True)"
```

Use `dry_run=true` to preview dispatcher selection without triggering children. To test the full
dispatch path without running Nextflow, use `child_dry_run=true` and an explicit
`max_jobs_per_kind`; the dispatcher triggers each child with `dry_run=true`, so it builds the
Nextflow command but does not execute it. In production, omit `child_dry_run` and it reads free
slots from each child deployment's work pool, using an idempotency key scoped to dispatcher run
and job ID.

TEST and PROD remain separate policy decisions. A successful TEST job does not automatically
promote to PROD, and source BAM/BAI or CRAM/CRAI files remain retained until accessions, receipts,
checksums, and archival verification have been recorded.

`nextflow_profile` is passed through to the pipeline. Production manifests should name the profile
that enables the pipeline's cluster executor (for example, `singularity,slurm`); the profile must
exist in the checked-out `ensembl-genes-nf` configuration.

Registry event publication, registry-driven readiness checks, QC gating, and handover dispatch are
deliberately outside this layer for now. The planner or a later integration layer must decide when
jobs are eligible and must not treat a successful local state file as a registry annotation event.
