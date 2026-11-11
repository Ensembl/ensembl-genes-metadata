# Prefect orchestration

`gb_prefect` runs Genebuild's Nextflow pipelines and Python scripts on Codon through
[Prefect](https://docs.prefect.io/): it schedules them, runs them as Slurm jobs, keeps their logs and
shows each run's status, logs and results in the Prefect UI.

These pages are for developers who know Python but not Prefect. Read this page first, then the
tutorial for what you want to add:

- [Add a Nextflow flow](add-nextflow-flow.md): run a Nextflow pipeline (this repo's or
  `ensembl-genes-nf`'s).
- [Add a script flow](add-script-flow.md): run a Python (or any) script, or Python code directly.

Then, as needed: [Deploying flows](deploying.md), [Credentials](credentials.md),
[Testing](testing.md), [Troubleshooting](troubleshooting.md), [Operations](operations.md) and
[Design notes](design-notes.md). The reference for every existing flow, its parameters and outputs
is the [`gb_prefect` README](https://github.com/Ensembl/ensembl-genes-metadata/blob/main/gb_prefect/README.md).

## Prefect in five minutes

| Prefect concept | What it is | Here |
|---|---|---|
| **Flow** | A Python function decorated with `@flow`: the unit you run, schedule and see in the UI. Its arguments are the run's parameters. | `gb_registry_flow`, `busco_gca_flow`, `gb_is_reference_flow` |
| **Task** | A function decorated with `@task`, called from a flow; each call is shown, logged and timed as a step of the flow run. | `run_nextflow_pipeline`, `run_script` |
| **Flow run** | One execution of a flow, with its parameters, state (Scheduled, Running, Completed, Failed, Crashed), logs and artifacts. | |
| **Deployment** | A flow registered on the Prefect server under a name, with default parameters, a work pool and job settings, so it can be started from the UI, the CLI, a schedule or another flow. | `gb_metadata_registry`, `busco-gca` |
| **Work pool** | Where a deployment's runs go to wait; its type decides *how* they run. | `genebuild-metadata` and `genebuild-busco` (`slurm-cli`: each run is a Slurm job), `genebuild-dispatch` (`process`: a plain process) |
| **Worker** | A long-running process that polls a work pool and starts its runs. It runs as a Unix user; everything it starts runs as that user. | `prefect worker start --pool genebuild-metadata --type slurm-cli` |
| **Secret block** | A credential stored on the server, loaded by name at run time, so it is never a deployment parameter. | `gb-metadata-db-params`, `gb-slack-params`, `nf-metro-token` |
| **Variable** | A plain value stored on the server. | `nf-metro-url` (turns live maps on) |
| **Artifact** | A summary attached to a flow run in the UI. | the run summary (command, status, log), the "Live pipeline map" link |

The worker clones the deployment's Git branch from GitHub for each run: **a deployment runs the code
pushed to its branch**, not your working copy.

## What happens when a flow runs

```text
 Prefect UI / CLI / schedule / another flow
        │  "run deployment gb_metadata_registry"
        ▼
 Prefect server (genebuild-metadata-prefect VM) ── the run waits in its work pool
        │
        ▼
 Worker for that pool (runs as a Unix user, e.g. genebuild)
        │  slurm-cli pool: writes <outdir>/slurm_jobs/slurm_cli_<run>.sh and sbatch-es it
        ▼
 Slurm job  ── setup_commands: module load nextflow; source <venv>/bin/activate
        │  clones the branch, runs the flow function (Python) in this job
        ▼
 Flow  ── resolves inputs (ENSCODE, credentials, folders), calls a task
        │
        ├─► run_nextflow_pipeline: nextflow run <ENSCODE>/<pipeline>/main.nf ...
        │       └─► Nextflow submits each pipeline process as its own Slurm job
        │
        └─► run_script: python <ENSCODE>/.../script.py ...   (inside this same job)

 Output: logs and records in the run's folder; logs, state and artifacts in the Prefect UI;
         optionally a live map of the Nextflow run (nf-metro, see Operations).
```

Two things follow from this picture, and most first-run problems come from them:

- **The worker's Unix user does everything**: it must be able to write the output directories, and
  the Slurm job, Nextflow and your script all run as that user.
- **A `slurm-cli` flow run *is* a Slurm job**: its `cpu`, `memory` and `time_limit` are that job's,
  and for a Nextflow flow that job is the Nextflow head job, which must outlive the whole pipeline.
  A `process` pool runs the flow on the worker's own machine, which suits only quick orchestration.

## The three kinds of flow

| Kind | Use it to | Built from | Example |
|---|---|---|---|
| **Nextflow flow** | run a Nextflow pipeline | a pipeline *spec* (`gb_prefect/pipelines/`) + `run_nextflow_pipeline` | `gb_registry_by_gca_flow` |
| **Script flow** | run a script or any command | a `ScriptRun` + `run_script` | `gb_is_reference_flow` |
| **Python flow** | do quick Python work in-process (DB queries, API calls, triggering other deployments) | your own `@task`s | the BUSCO dispatchers |

## Where the code lives

| Path in `gb_prefect/` | What goes there |
|---|---|
| `flows/<group>/` | Flows, grouped by source (`assembly_metadata/`, `statistics/`). |
| `pipelines/` | One spec module per Nextflow pipeline (mode), e.g. `assembly_metadata.py`, `statistics_busco_genome.py`. Plain Python. |
| `tasks/` | Tasks: the generic `nextflow.py` and `script.py`, plus flow-specific steps. |
| `models/` | Data classes: `NextflowRun`, `ScriptRun`, results and errors, `PipelineCredentials`. |
| `nextflow/` | Building and running commands, nf-metro registration. Plain Python. |
| `utils/` | Helpers: ENSCODE, credentials, folders, input files, artifacts, deployment. |
| `setup/` | `deploy_flows.py` (all deployments) and `create_secrets.py` (Secret blocks). |
| `worker/` | The `slurm-cli` worker type. |
