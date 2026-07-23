# Genebuild Metadata Prefect

`gb_prefect` is a Python package of [Prefect](https://www.prefect.io/) flows and tasks that submit Genebuild's Nextflow pipelines (and standalone scripts) as SLURM jobs, capture their logs, and optionally publish a run summary as a Prefect artifact.

Each flow can be run two ways:

- **As a standalone script** (`python gb_prefect/flows/<flow>.py --...`) — no Prefect server required.
- **As a Prefect flow**, imported and called from Python (e.g. from another flow, a Prefect deployment, or an interactive session).

---

## Table of Contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Package Layout](#package-layout)
- [Flows](#flows)
- [Common Behavior](#common-behavior)
- [Outputs](#outputs)

---

## Requirements

- Python 3.9+ and [Prefect](https://www.prefect.io/)
- Access to a SLURM cluster (flows submit jobs via `sbatch --wait`)
- The `nextflow/24.10.3` environment module available on the cluster (for all flows except `gb_is_reference`)
- An `ENSCODE` directory containing the Ensembl repositories referenced by each flow (see [Flows](#flows) below) — pass it explicitly via `--enscode`/`enscode=` or set the `ENSCODE` environment variable
- A Python virtual environment (`asm_venv`) on the cluster with the dependencies needed by the pipeline/script being run.  

## Installation

From the repository root:

```bash
pip install .
```

This installs `gb_prefect` as part of the `ensembl-genes-metadata` package (see [pyproject.toml](../pyproject.toml)).

## Package Layout

| Path | Contents |
|------|----------|
| [flows/](flows) | Prefect `@flow`-decorated entry points, each with a CLI (`argparse`) for standalone use |
| [tasks/](tasks) | Prefect `@task`-decorated functions that build a SLURM/sbatch script, submit it, and capture the result |
| [utils/](utils) | Shared helpers: ENSCODE resolution, logging, shell execution, artifact publishing, CSV splitting |
| [models/](models) | Reserved for shared data models (currently empty) |

---

## Flows

| Flow | Module | Runs | Purpose |
|------|--------|------|---------|
| `gb_registry_flow` | [flows/gb_registry.py](flows/gb_registry.py) | [pipelines/assembly_metadata](../pipelines/assembly_metadata) | Register new assemblies released on a given date |
| `gb_registry_update_flow` | [flows/gb_registry_update.py](flows/gb_registry_update.py) | [pipelines/assembly_metadata_update](../pipelines/assembly_metadata_update) | Check/update metadata for a list of existing GCA accessions |
| `genome_busco_flow` | [flows/gb_busco_genome.py](flows/gb_busco_genome.py) | `ensembl-genes-nf/pipelines/statistics` | Run BUSCO genome statistics for a single CSV file |
| `genome_busco_master_flow` | [flows/gb_busco_genome_single_bulk.py](flows/gb_busco_genome_single_bulk.py) | `ensembl-genes-nf/pipelines/statistics` | Split a CSV into one row per file and run `genome_busco_flow`'s task in parallel for each |
| `gb_is_reference_flow` | [flows/gb_is_reference.py](flows/gb_is_reference.py) | [src/python/is_reference.py](../src/python/is_reference.py) | Check whether GCA accessions are the NCBI reference assembly for their taxon |

### `gb_registry_flow`

```bash
python gb_prefect/flows/gb_registry.py \
    --date 01-15-2026 \
    --outdir /path/to/output \
    --enscode $ENSCODE \
    --asm_venv /path/to/venv \
    [--dry-run]
```

| Argument | Required | Description |
|----------|----------|-------------|
| `--date` | Yes | Date to register assemblies for. Format: `MM-DD-YYYY` |
| `--outdir` | Yes | Base output directory (a `<date>` subdirectory is created under it) |
| `--enscode` | Yes | Path to the `ENSCODE` directory |
| `--asm_venv` | Yes | Path to the virtual environment used to run the Nextflow pipeline |
| `--dry-run` | No | Build the sbatch script and log it, but don't submit the job |

### `gb_registry_update_flow`

```bash
python gb_prefect/flows/gb_registry_update.py \
    --gca-list /path/to/gca_list.txt \
    --outdir /path/to/output \
    --asm_venv /path/to/venv \
    --enscode $ENSCODE \
    [--date YYYY-MM-DD] \
    [--slack-report] \
    [--dry-run]
```

| Argument | Required | Description |
|----------|----------|-------------|
| `--gca-list` | Yes | Path to a file listing the GCA accessions to check/update |
| `--outdir` | Yes | Base output directory (an `asm_update_<date>` subdirectory is created under it) |
| `--asm_venv` | Yes | Path to the virtual environment used to run the Nextflow pipeline |
| `--enscode` | Yes | Path to the `ENSCODE` directory |
| `--date` | No | Defaults to today (`YYYY-MM-DD`) if not set |
| `--slack-report` | No | Enables Slack reporting in the underlying pipeline |
| `--dry-run` | No | Build the sbatch script and log it, but don't submit the job |

> **Note:** the sbatch script this task generates currently passes a hardcoded `--screen_date 2026-01-15` to the pipeline alongside `--gca_input`/`--gca_list`. Since [assembly_metadata_update](../pipelines/assembly_metadata_update) treats `--screen_date` and `--gca_input` as mutually exclusive, this value is unused when running against a GCA list — but it's worth updating or removing if this flow is reused after that date.

### `genome_busco_flow`

```bash
python gb_prefect/flows/gb_busco_genome.py \
    --csv-file /path/to/genomes.csv \
    --outdir /path/to/output \
    --asm-venv /path/to/venv \
    [--enscode $ENSCODE] \
    [--dry-run] \
    [--create-artifact]
```

| Argument | Required | Description |
|----------|----------|-------------|
| `--csv-file` | Yes | CSV file of genomes to run BUSCO statistics on |
| `--outdir` | Yes | Output directory |
| `--asm-venv` | Yes | Path to the virtual environment used to run the Nextflow pipeline |
| `--enscode` | No | Path to the `ENSCODE` directory; falls back to the `ENSCODE` environment variable |
| `--dry-run` | No | Build the sbatch script and log it, but don't submit the job |
| `--create-artifact` | No | Publish a Prefect markdown artifact summarizing the run |

### `genome_busco_master_flow`

Same arguments as `genome_busco_flow` (see [flows/gb_busco_genome_single_bulk.py](flows/gb_busco_genome_single_bulk.py)), except `--csv-file` is split into one row per file first, and `genome_busco_flow`'s underlying task is submitted once per row **in parallel**, returning a list of results.

### `gb_is_reference_flow`

```bash
python gb_prefect/flows/gb_is_reference.py \
    --file-path /path/to/gca_list.txt \
    --output-path /path/to/output \
    --enscode $ENSCODE \
    --asm_venv /path/to/venv
```

| Argument | Required | Description |
|----------|----------|-------------|
| `--file-path` | Yes | File of GCA accessions to check (passed through to [is_reference.py](../src/python/is_reference.py)) |
| `--output-path` | Yes | Directory where the output CSV and logs will be saved |
| `--enscode` | Yes | Path to the `ENSCODE` directory |
| `--asm_venv` | Yes | Path to the virtual environment used to run `is_reference.py` |

Unlike the other flows, this one has no `--dry-run` or `--create-artifact` option, and raises `RuntimeError` (rather than returning a Prefect `Failed` state) if the SLURM job fails.

---

## Common Behavior

- **ENSCODE resolution** ([utils/enscode_utils.py](utils/enscode_utils.py)): for every flow except `gb_is_reference`, `enscode` is resolved from the argument, then the `ENSCODE` environment variable, then — only when `dry_run=True` — a `<ENSCODE>` placeholder so the generated command can still be previewed. If none is available and `dry_run=False`, a `ValueError` is raised.
- **Dry runs**: passing `dry_run=True` (or `--dry-run`) builds and logs the sbatch script and command file without calling `sbatch`. `returncode` is `0` and no SLURM job ID is produced.
- **Logging** ([utils/logging_utils.py](utils/logging_utils.py)): every task appends timestamped `INFO` lines to a `log_flow_*.log` file in the output directory, including the full sbatch script and, once the job finishes, the contents of SLURM's own `.out`/`.err` files.

## Outputs

For each run, in the given output directory:

| File | Produced by | Description |
|------|-------------|-------------|
| `log_flow_*.log` | every task | Timestamped log of the run, including the sbatch script and SLURM output |
| `*_command_*.sh` | every task | The generated, executable sbatch script that was (or would have been) submitted |
| `slurm_<job_id>.out` / `.err` (naming varies per task) | SLURM | Raw stdout/stderr from the submitted job |
| Prefect markdown artifact | `register_assemblies`, `update_assemblies`, `run_nextflow_busco` (when `create_artifact=True`, the default) | Run summary (status, command, log excerpt) visible in the Prefect UI |

Return values differ slightly by task:

- `register_assemblies` / `update_assemblies` return a Prefect `Completed`/`Failed` state wrapping a result dict (`returncode`, `command`, `command_file`, `log_file`, `slurm_job_id`, `pipeline_run_date`, `pipeline_ran`, `dry_run`).
- `run_nextflow_busco` returns the same kind of result dict directly, without wrapping it in a Prefect state.
- `is_reference` returns the output CSV path as a string on success, or raises `RuntimeError` on failure.
