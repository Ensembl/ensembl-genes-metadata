# Genebuild Metadata Prefect

`gb_prefect` is a Python package of [Prefect](https://www.prefect.io/) flows and tasks that run Genebuild's Nextflow pipelines (and standalone scripts) on SLURM, capture their logs, and publish a run summary as a Prefect artifact.

Nextflow flows run on the [`slurm-cli` worker](#slurm-cli-worker): each flow run is one Slurm job, and inside it the generic task [`run_nextflow_pipeline`](tasks/nextflow.py) launches Nextflow directly. A flow only describes *what* to run, through a pipeline spec in [pipelines/](pipelines); see [Nextflow Flows](#nextflow-flows).

Flows are run through their Prefect deployments (see [Running a Flow](#running-a-flow)); there is no per-flow command line.

**New to Prefect, or adding a flow?** Start with the guides on the docs site: [Prefect orchestration](https://ensembl.github.io/ensembl-genes-metadata/prefect/) (concepts and architecture), [Add a Nextflow flow](https://ensembl.github.io/ensembl-genes-metadata/prefect/add-nextflow-flow/), [Add a script flow](https://ensembl.github.io/ensembl-genes-metadata/prefect/add-script-flow/), [Testing](https://ensembl.github.io/ensembl-genes-metadata/prefect/testing/) and [Troubleshooting](https://ensembl.github.io/ensembl-genes-metadata/prefect/troubleshooting/) (sources in [`docs/prefect/`](../docs/prefect)). This README is the reference for the existing flows.

---

## Table of Contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Package Layout](#package-layout)
- [Flows](#flows)
- [Setting up the Prefect Server](#setting-up-the-prefect-server)
- [Running a Flow](#running-a-flow)
- [Nextflow Flows](#nextflow-flows)
- [Script Flows](#script-flows)
- [Live Maps](#live-maps)
- [Common Behavior](#common-behavior)
- [Outputs](#outputs)
- [`slurm-cli` Worker](#slurm-cli-worker)

---

## Requirements

- Python 3.9+ and [Prefect](https://www.prefect.io/)
- Access to a SLURM cluster
- A `nextflow` environment module on the cluster, loaded by the work pool's `setup_commands` (see [`deploy_flows.py`](setup/deploy_flows.py))
- An `ENSCODE` directory containing the Ensembl repositories referenced by each flow (see [Flows](#flows) below) — pass it as the `enscode` parameter or set the `ENSCODE` environment variable
- A Python virtual environment (`asm_venv`) on the cluster with this package, `prefect` and the pipelines' dependencies installed, activated by the work pool's `setup_commands`

## Installation

From the repository root:

```bash
pip install .
```

This installs `gb_prefect` as part of the `ensembl-genes-metadata` package (see [pyproject.toml](../pyproject.toml)).

## Package Layout

| Path | Contents |
|------|----------|
| [flows/](flows) | Prefect `@flow`-decorated entry points, grouped by source pipeline ([`assembly_metadata/`](flows/assembly_metadata), [`statistics/`](flows/statistics)); their parameters are documented in each flow's docstring |
| [pipelines/](pipelines) | Pipeline specs, one module per Nextflow pipeline (mode): build a `NextflowRun` from the flow's inputs; no Prefect |
| [nextflow/](nextflow) | Turning a `NextflowRun` into its command, params file and input CSV, and running it with live output; no Prefect |
| [tasks/](tasks) | Prefect `@task`s: the generic [`run_nextflow_pipeline`](tasks/nextflow.py) and [`run_script`](tasks/script.py), post-run steps ([`statistics_cleanup`](tasks/statistics_cleanup.py)), and BUSCO dispatch |
| [models/](models) | Shared data models: `NextflowRun` / `NextflowResult` / `NextflowRunError`, `ScriptRun` / `ScriptResult` / `ScriptRunError`, `PipelineCredentials`, `TaskRunOptions` |
| [utils/](utils) | Shared helpers: ENSCODE, credentials, output folders and input files, logging, artifacts, CSV |
| [setup/](setup) | One-time Prefect server setup scripts: [`create_secrets.py`](setup/create_secrets.py) (Secret blocks) and [`deploy_flows.py`](setup/deploy_flows.py) (deployments); see [Setting up the Prefect Server](#setting-up-the-prefect-server) |
| [worker/](worker) | The `slurm-cli` Prefect worker type |

---

## Flows

| Flow | Module | Runs | Purpose |
|------|--------|------|---------|
| `gb_registry_flow` | [flows/assembly_metadata/gb_registry.py](flows/assembly_metadata/gb_registry.py) | [pipelines/assembly_metadata](../pipelines/assembly_metadata) | Register new assemblies released after a given date (default: 60 days ago) |
| `gb_registry_by_gca_flow` | [flows/assembly_metadata/gb_registry_by_gca.py](flows/assembly_metadata/gb_registry_by_gca.py) | [pipelines/assembly_metadata](../pipelines/assembly_metadata) | Register the GCA accessions listed in a file |
| `gb_registry_update_flow` | [flows/assembly_metadata/gb_registry_update.py](flows/assembly_metadata/gb_registry_update.py) | [pipelines/assembly_metadata_update](../pipelines/assembly_metadata_update) | Check/update metadata for a list of existing GCA accessions |
| `gb_registry_update_by_date_flow` | [flows/assembly_metadata/gb_registry_update_by_date.py](flows/assembly_metadata/gb_registry_update_by_date.py) | [pipelines/assembly_metadata_update](../pipelines/assembly_metadata_update) | Check/update metadata for high-priority assemblies released after a date |
| `busco_gca_flow` | [flows/statistics/gb_busco_gca.py](flows/statistics/gb_busco_gca.py) | `ensembl-genes-nf/pipelines/statistics` | Run BUSCO genome statistics for one GCA (triggered per GCA by the BUSCO dispatchers) |
| `busco_dispatch_flow` | [flows/statistics/gb_busco_dispatch.py](flows/statistics/gb_busco_dispatch.py) | `busco_gca_flow` | Trigger `busco_gca_flow` for each GCA of a CSV (`gca,taxon_id,busco_dataset`) |
| `busco_dispatch_from_list_flow` | [flows/statistics/gb_busco_dispatch_from_list.py](flows/statistics/gb_busco_dispatch_from_list.py) | `busco_gca_flow` | Trigger `busco_gca_flow` for GCAs given as a list or a file |
| `busco_dispatch_automatic_flow` | [flows/statistics/gb_busco_dispatch_automatic.py](flows/statistics/gb_busco_dispatch_automatic.py) | `busco_gca_flow` | Pick pending BUSCO candidates from the registry, sized to the BUSCO pool's free slots |
| `gb_is_reference_flow` | [flows/gb_is_reference.py](flows/gb_is_reference.py) | [src/python/is_reference.py](../src/python/is_reference.py) | Check whether GCA accessions are the NCBI reference assembly for their taxon |

The registry and update flows and `busco_gca_flow` are [Nextflow flows](#nextflow-flows) and `gb_is_reference_flow` is a [script flow](#script-flows); they run on a `slurm-cli` work pool (`genebuild-metadata` and `genebuild-busco`). The BUSCO dispatchers run no Nextflow and use the `genebuild-dispatch` process pool. Every flow's parameters, and their defaults, are in its docstring; the main ones are below.

### `gb_registry_flow`

Screens NCBI for assemblies released after a date and registers the ones missing from the metadata DB.

| Parameter | Description |
|-----------|-------------|
| `outdir` | Base output directory; the run goes to `<outdir>/<folder_name>` |
| `date` | Screen for assemblies released after this date (`YYYY-MM-DD`, converted to the pipeline's `MM/DD/YYYY`). Default: 60 days before today |
| `folder_name` | Output folder under `outdir` (a single name, not a path). Default: `asm_registry_<today>` |
| `enscode` | `ENSCODE` directory. Default: the `ENSCODE` environment variable |
| `dry_run` | Write and log the command file without running Nextflow or loading credentials |
| `credentials`, `metadata_secret_block` | DB credentials: explicit `PipelineCredentials`, or loaded from the Secret block (default `gb-metadata-db-params`) |

### `gb_registry_by_gca_flow`

Registers the GCA accessions listed in `gca_list` (`--add_gca` mode, the whole file as given), without screening NCBI. `gca_list` must be an absolute path to an existing file (existence is not checked on a dry run). `folder_name` defaults to `asm_registry_gca_<today>`; the other parameters are as for `gb_registry_flow`.

### `gb_registry_update_flow`

Checks/updates the metadata of the GCA accessions listed in `gca_list` (an absolute path; the whole file, as given). `folder_name` defaults to `asm_update_<today>`. `slack_report` (default true) sends a Slack report at the end, with Slack credentials from `credentials` or the `gb-slack-params` Secret block (`slack_secret_block`). The other parameters are as for `gb_registry_flow`.

### `gb_registry_update_by_date_flow`

Screens the metadata DB for high-priority assemblies released after `screen_date` (`YYYY-MM-DD`) and checks/updates them, instead of using a GCA list. `folder_name` defaults to `asm_update_screen_<screen_date>_<today>`; the other parameters are as for `gb_registry_update_flow`.

### `busco_gca_flow`

Runs the statistics pipeline's BUSCO genome mode for one `gca` / `taxon_id` under `<outdir>/<gca>`; `busco_dataset` picks the lineage (default: chosen by the pipeline). `update_registry` loads the results into the metadata DB and marks the GCA done; `post_run_cleanup` (set by the automatic dispatcher) removes the run folder on success, or marks the GCA failed and keeps the logs on failure.

### `gb_is_reference_flow`

Checks whether the GCA accessions in `file_path` (an absolute path; one accession per line) are the NCBI reference assembly for their taxon, with [is_reference.py](../src/python/is_reference.py) from ENSCODE. The run goes to `<outdir>/<folder_name>` (`folder_name` defaults to `is_reference_<today>`), and the flow returns the path of the CSV it writes there, `<file_path stem>_output.csv`. It needs no credentials (NCBI's public API); `enscode` and `dry_run` are as for `gb_registry_flow`.

---

## Setting up the Prefect Server

One-time steps, in this order, with `PREFECT_API_URL` pointing at the Prefect server and the venv that has this package installed (`pip install .`) active. Every step can be re-run.

1. **Create the work pools.** The `slurm-cli` type comes from this package (see [`slurm-cli` Worker](#slurm-cli-worker)):

   ```bash
   prefect work-pool create genebuild-metadata --type slurm-cli
   ```
   ```bash
   prefect work-pool create genebuild-busco --type slurm-cli
   ```
   ```bash
   prefect work-pool update genebuild-busco --concurrency-limit 10
   ```
   ```bash
   prefect work-pool create genebuild-dispatch --type process
   ```

   `genebuild-metadata` runs the registry/update flows (no pool limit; each deployment allows one run at a time). `genebuild-busco` runs `busco_gca_flow`; its concurrency limit is how many BUSCO runs go at once, and the automatic dispatcher refuses to run without one. `genebuild-dispatch` is a process pool for the BUSCO dispatchers, which run no Nextflow.

2. **Start a worker for each pool**, as the user that owns the output directories: a `slurm-cli` worker writes each run's sbatch script there, and the Slurm job, Nextflow included, runs as that user.

   ```bash
   prefect worker start --pool genebuild-metadata --type slurm-cli
   ```
   ```bash
   prefect worker start --pool genebuild-dispatch --type process
   ```

   (and the same for `genebuild-busco`).

3. **Create the Secret blocks** with the credentials, from JSON files (see [`create_secrets.py`](setup/create_secrets.py) for their format):

   ```bash
   python gb_prefect/setup/create_secrets.py --metadata-json ~/secrets/metadata_params.json --slack-json ~/secrets/slack_params.json
   ```

   `gb-metadata-db-params` is used by every flow that reads or writes the metadata DB; `gb-slack-params` only by the update flows' Slack report. Add `--overwrite` to replace existing blocks.

4. **Deploy the flows** with [`deploy_flows.py`](setup/deploy_flows.py) (commit and push the branch first: the worker clones it at run time). `--dry-run` shows what would be deployed:

   ```bash
   python gb_prefect/setup/deploy_flows.py --branch main --outdir /path/to/asm_output --flow-outdir busco_gca=/path/to/busco_output --enscode /path/to/enscode --asm-venv /path/to/asm_venv
   ```

   Schedules (e.g. an hourly `busco-dispatch-automatic`) are set in the Prefect UI; check they are still there after redeploying.

5. **Optional: live maps.** Run an nf-metro server where both the Codon compute nodes and your browser can reach it (we use the Prefect VM, port 8081), with a random token (`python3 -c "import secrets; print(secrets.token_urlsafe(32))"`):

   ```bash
   nf-metro serve-multi --host 0.0.0.0 --port 8081 --token <token>
   ```

   Then point Prefect at it, and store the same token (prompted, not echoed):

   ```bash
   prefect variable set nf-metro-url http://genebuild-metadata-prefect.ebi.ac.uk:8081
   ```
   ```bash
   python gb_prefect/setup/create_secrets.py --metro-token
   ```

   See [Live Maps](#live-maps).

## Running a Flow

With the server set up, start runs from the Prefect UI or the CLI, overriding any parameter:

```bash
prefect deployment run 'gb_registry/gb_metadata_registry' --param dry_run=true
```

```bash
prefect deployment run 'gb_registry_update_by_date/gb_metadata_update_per_date' --param screen_date=2026-01-15
```

A flow can also be called directly from Python, e.g. for a dry run that shows the generated Nextflow command without a deployment (no Secret block is needed on a dry run):

```bash
python -c "from gb_prefect.flows.assembly_metadata.gb_registry import gb_registry_flow as f; f(outdir='/path/to/output', dry_run=True)"
```

Called directly, a Nextflow flow launches Nextflow in the current shell, so for a real run do it from a Slurm job or interactive session with `nextflow` and the venv loaded. Pass credentials as `credentials=PipelineCredentials(metadata_params_string='{...}')` (from `gb_prefect.models.pipeline_options`) to avoid the Secret blocks.

---

## Nextflow Flows

A Nextflow flow builds a [`NextflowRun`](models/nextflow.py) with a pipeline spec from [pipelines/](pipelines) and passes it to the generic task [`run_nextflow_pipeline`](tasks/nextflow.py), which:

1. writes the input CSV (if the spec builds one) and an owner-only (`0600`) `-params-file` with the credentials, so they never appear on the command line, in the logs or in the artifact;
2. runs `nextflow` directly (no shell) from the run's folder, streaming its output line by line to the log file and the Prefect logs;
3. writes the equivalent command as a record file you can re-run by hand (`bash nextflow_<name>_command.sh`), and publishes a Prefect artifact;
4. raises `NextflowRunError` (carrying a `NextflowResult`) if Nextflow exits non-zero, **or exits 0 but its output shows tasks that failed and were ignored** (a pipeline with `errorStrategy 'ignore'` reports "Successful" although those tasks produced nothing; set `NextflowRun.fail_on_ignored_errors=False` to accept them); otherwise returns the `NextflowResult`.

Terminal colour codes are stripped from Nextflow's output, and the artifact shows only the current run's part of the log (the log file itself is appended to across reruns in the same folder).

When the flow run is a Prefect retry, `-resume` is added, so Nextflow reuses the cached tasks in the same folder. DAG files are timestamped (`*_dag_<YYYYmmdd-HHMMSS>.mmd`) so reruns in the same folder don't clash.

Adding a pipeline means adding a spec module in [pipelines/](pipelines) and a small flow that resolves its inputs (ENSCODE, credentials, folder) and calls the task.

## Script Flows

A script flow runs one command, typically a Python script from ENSCODE, with the generic task [`run_script`](tasks/script.py), the counterpart of `run_nextflow_pipeline`. The flow builds a [`ScriptRun`](models/script.py):

```python
run = ScriptRun(
    name="is_reference",                      # names the log, command record and artifact
    argv=[sys.executable, f"{enscode}/ensembl-genes-metadata/src/python/is_reference.py",
          "--file", file_path, "--output", output_file],
    outdir=run_outdir,                        # the command runs from here
    env={},                                   # extra environment variables
    secret_env={},                            # credentials: redacted, never in the record
)
run_script(run, TaskRunOptions(dry_run=dry_run))
```

`run_script` runs the command inside the flow run's own Slurm job (deploy the flow to a `slurm-cli` pool and give it the job's resources in [`deploy_flows.py`](setup/deploy_flows.py)); it does not submit a job of its own. `sys.executable` is the venv the pool's `setup_commands` activate. As with Nextflow, the command runs directly (no shell), its output is streamed to `<outdir>/script_<name>.log` and the Prefect logs, an equivalent bash script is written to `script_<name>_command.sh` as a record, an artifact shows this run's part of the log, and a non-zero exit raises `ScriptRunError` (carrying a `ScriptResult`). Pass credentials in `secret_env`, never in `argv`.

For quick Python work (DB queries, calling the API), a flow can also call functions directly, as the BUSCO dispatchers do with `gb_metadata`; such flows can use the `genebuild-dispatch` process pool (`runs_nextflow=False` in `deploy_flows.py`) instead of taking a Slurm slot.

## Live Maps

When the Prefect Variable `nf-metro-url` is set, every real run of a Nextflow flow whose pipeline has an nf-metro map (`NextflowRun.metro_map`) gets a live map of its progress. Maps live in the pipeline's `assets/` folder: `metro_map.mmd` for this repo's pipelines, and one map per mode for a pipeline with several (`metro_map_<mode>.mmd`, e.g. `ensembl-genes-nf/pipelines/statistics/assets/metro_map_busco_genome.mmd`):

1. the task registers the run with the [nf-metro](https://github.com/seqeralabs/nf-metro) server ([`nextflow/metro.py`](nextflow/metro.py)), sending the map from ENSCODE;
2. Nextflow is launched with `-with-weblog`, so its events reach the server and light up each process's station;
3. a **"Live pipeline map"** link artifact on the flow run opens the map.

Nothing here can fail a run: if the server is down, the token is wrong or the map is missing, the run logs a warning and goes on without a live map. Dry runs register nothing. Delete the Variable to switch live maps off. The token (Secret block `nf-metro-token`) is redacted from the command record, log and artifacts; Nextflow itself still writes its full command line, token included, to `.nextflow.log`, so use a token that protects nothing else.

The server keeps runs in memory only (its newest 100): a link stops working after a server restart. The map is read from ENSCODE, so the ENSCODE checkout must contain the maps (`git pull` it). Edit maps and re-render the diagrams with [`src/bash/render_metro_maps.sh`](../src/bash/render_metro_maps.sh).

## Common Behavior

- **ENSCODE resolution** ([utils/enscode_utils.py](utils/enscode_utils.py)): for every flow except `gb_is_reference`, `enscode` is resolved from the argument, then the `ENSCODE` environment variable, then — only when `dry_run=True` — a `<ENSCODE>` placeholder so the generated command can still be previewed. If none is available and `dry_run=False`, a `ValueError` is raised.
- **Credentials** ([utils/credentials_utils.py](utils/credentials_utils.py)): passed explicitly (`credentials=PipelineCredentials(...)`) or loaded from Prefect Secret blocks (the default, for deployments). Not loaded on a dry run. In a flow:

  ```python
  from gb_prefect.utils.credentials_utils import PipelineCredentials, resolve_credentials

  credentials = resolve_credentials(credentials, metadata_secret_block=metadata_secret_block)
  credentials.metadata_params_string   # JSON string, passed to a pipeline
  credentials.metadata_db_params       # dict, for code connecting to the DB itself
  ```

  Add `slack_report=True` to also load the Slack credentials, and `dry_run=dry_run` to load nothing on a dry run. A flow's `metadata_secret_block` / `slack_secret_block` parameters default to `None`, meaning the standard blocks (`gb-metadata-db-params`, `gb-slack-params`).
- **Dry runs**: passing `dry_run=True` writes and logs the command file without running anything. `returncode` is `0`.

## Outputs

For each Nextflow flow run, in its output folder:

| File | Description |
|------|-------------|
| `nextflow_<name>.log` | Timestamped log of the run: the command, then Nextflow's full output |
| `nextflow_<name>_command.sh` | Record of the command that was (or would have been) run; re-runnable by hand |
| `nextflow_<name>_params.json` | Credentials passed with `-params-file`, owner-only (not written on a dry run) |
| `<name>.csv` | Input CSV, when the spec builds one (e.g. one row per GCA for BUSCO) |
| `*_dag_<timestamp>.mmd`, `report-<timestamp>.html` | Nextflow DAG (Mermaid) and execution report, when the spec asks for them |
| Prefect markdown artifact | Run summary (status, command, last 300 log lines) in the Prefect UI |

Each deployment's Slurm job also writes, in `<outdir>/slurm_jobs/` (`--job-dir` when deploying): the worker's sbatch script `slurm_cli_<flow run id>.sh` and Slurm's `slurm_<job id>.out` / `.err`. Look there first when a run crashes before the flow starts.

A script flow writes `script_<name>.log` and `script_<name>_command.sh` instead, plus whatever the script itself writes.

---

## `slurm-cli` Worker

`gb_prefect.worker.SlurmCliWorker` is a custom Prefect worker type (`type: slurm-cli`) that submits an entire flow run as a single Slurm job via the `sbatch` CLI. It's registered under `[project.entry-points."prefect.collections"]` in [pyproject.toml](../pyproject.toml), so once the package is installed it's available to `prefect work-pool create` / `prefect worker start` like any built-in worker type.

As with any Prefect worker, the **whole flow run** — including any task that shells out to Nextflow — executes inside one Slurm job; tasks are not separate Slurm jobs.

### Job configuration fields

| Field | Description |
|---|---|
| `cpu` | CPU count (`--cpus-per-task`) |
| `memory` | Memory in GB (`--mem`) |
| `partition` | Slurm partition (`--partition`), optional |
| `account` | Slurm account (`--account`), optional |
| `qos` | Slurm QOS (`--qos`), optional |
| `time_limit` | Wall time in hours (`--time`) |
| `working_dir` | Directory the job runs from; also where the generated sbatch script and `slurm_%j.out`/`.err` files are written. [`deploy_flows.py`](setup/deploy_flows.py) sets it to `<flow outdir>/slurm_jobs` (`--job-dir` to rename) |
| `setup_commands` | Shell commands run before the flow-run command, e.g. `module load nextflow/25.04.6`, `source ~/asm_venv/bin/activate` |

**Important**: the environment produced by `setup_commands` must have this package (`gb_prefect`) and `prefect` importable — the flow run's own Python process executes inside the submitted Slurm job, not just any Nextflow subprocess it launches.

### Creating a pool

See [Setting up the Prefect Server](#setting-up-the-prefect-server) for the pools the flows use. Job configuration values (`cpu`, `memory`, `setup_commands`, etc.) are then set on the pool's base job template or overridden per-deployment, the same as any other Prefect worker type.

### Current scope

This is an initial implementation. `sbatch` submission is retried with backoff (`GB_PREFECT_WORKER_MAX_ATTEMPTS`, default 3, plus `GB_PREFECT_WORKER_RETRY_MIN_DELAY_SECONDS`/`..._JITTER_SECONDS` env vars) for transient `slurmctld` errors.
