# Deploying flows

A flow runs on the server only through a **deployment**: the flow registered under a name, with a
work pool, default parameters and job settings. All deployments are described in one file,
[`gb_prefect/setup/deploy_flows.py`](https://github.com/Ensembl/ensembl-genes-metadata/blob/main/gb_prefect/setup/deploy_flows.py),
and created or updated by running it.

## Choosing a work pool

| Pool | Type | Runs | Use it for |
|---|---|---|---|
| `genebuild-metadata` | `slurm-cli` | each flow run as a Slurm job | registry/update flows, script flows |
| `genebuild-busco` | `slurm-cli` | each flow run as a Slurm job; the pool's concurrency limit caps how many run at once | `busco_gca` |
| `genebuild-dispatch` | `process` | each flow run as a process on the worker's machine | quick orchestration (the BUSCO dispatchers) |

- Anything that runs Nextflow or does real work goes on a **`slurm-cli`** pool. A Nextflow flow on a
  `process` pool runs the Nextflow head job on the worker's machine, where it can be killed
  (`return code -9`, see [Troubleshooting](troubleshooting.md)).
- A **`process`** pool suits flows that only make a few calls and must not take a Slurm slot.
- A new kind of work with its own throughput (e.g. OMArk) gets **its own pool**, so its concurrency
  limit doesn't compete with others. Add it to `POOLS` (below) and create it on the server
  ([Operations](operations.md)).

## A `FlowDeployment`

Each entry of `FLOWS` in `deploy_flows.py` is one deployment:

```python
    "register": FlowDeployment(
        name="gb_metadata_registry",           # the deployment's name on the server
        entrypoint="gb_prefect/flows/assembly_metadata/gb_registry.py:gb_registry_flow",
        work_pool="metadata",                   # a key of POOLS, not the pool's name
        description="Assembly metadata register",
        tags=["metadata_registry"],
        cpu=1, memory=4, time_limit=2,          # the Slurm job: CPUs, GB, hours
        concurrency_limit=1,                    # at most one run of this deployment at a time
    ),
```

| Field | Default | Notes |
|---|---|---|
| `name` | | Deployment name; runs are addressed as `<flow name>/<deployment name>`, e.g. `gb_registry/gb_metadata_registry`. |
| `entrypoint` | | `<path from the repo root>:<flow function>`. |
| `work_pool` | | A key of `POOLS`, e.g. `"metadata"`. |
| `cpu`, `memory`, `time_limit` | 1, 4, 2 | The flow run's Slurm job (CPUs, GB, hours). For a Nextflow flow this is the **head job**, which must outlive the whole pipeline; the pipeline's processes are separate Slurm jobs sized by the pipeline. |
| `concurrency_limit` | none | Max simultaneous runs of this deployment; extra runs wait. Use `1` when two runs must never write to the same DB at once. |
| `parameters` | | Extra default parameters, on top of `outdir` and `enscode`. |
| `separate_outdir` | `False` | `True` when the flow's output normally lives apart from `--outdir` (BUSCO): a warning is printed if it isn't given its own. |
| `runs_nextflow` | `True` | `True`: a Slurm-job flow (Nextflow *or* script) that gets `outdir`/`enscode` parameters and the Slurm job settings. `False`: a process-pool flow that gets only its own `parameters`, and runs with `<asm-venv>/bin/python`. |

Never put credentials in `parameters`: deployment parameters are stored and shown in plain text.
Flows load credentials from Secret blocks ([Credentials](credentials.md)).

`POOLS` maps each key to a pool name. Each key becomes a command-line option (`--metadata-pool`,
`--busco-pool`, `--dispatch-pool`) to deploy to a different pool, e.g. your own test pool.

## Running the deploy script

From the repository root, with the venv that has this package installed and `PREFECT_API_URL`
pointing at the server. **Commit and push the branch first**: the worker clones it for every run.

```bash
python gb_prefect/setup/deploy_flows.py --branch main --outdir /path/to/asm_output --flow-outdir busco_gca=/path/to/busco_output --enscode /path/to/enscode --asm-venv /path/to/asm_venv
```

| Option | Notes |
|---|---|
| `--branch` | The branch the worker clones (default `main`). |
| `--flows KEY ...` | Deploy only these (default: all). |
| `--name-suffix _test` | Deploy copies (`gb_metadata_registry_test`) instead of replacing the production deployments. |
| `--outdir` | Base output directory of every Slurm-job flow; also where each job's files go (below). |
| `--flow-outdir KEY=PATH` | A different output directory for one flow (repeatable), e.g. `busco_gca=/path/busco`. |
| `--job-dir NAME` | The folder for each job's sbatch script and Slurm `.out`/`.err` (default `slurm_jobs`). |
| `--enscode`, `--asm-venv`, `--nextflow-module` | What the job uses: ENSCODE, the venv it activates, the Nextflow module it loads. |
| `--metadata-pool` / `--busco-pool` / `--dispatch-pool POOL` | Deploy that pool's flows to another pool. |
| `--work-pool POOL` | Deploy every selected Slurm-job flow to one pool (process-pool flows keep theirs). |
| `--cpu`, `--memory`, `--time-limit`, `--partition` | Override the Slurm settings of every selected flow. |
| `--dry-run` | Print what would be deployed; deploy nothing. Needs no server. |

The worker writes each run's sbatch script (`slurm_cli_<flow run id>.sh`) and Slurm's
`slurm_<job id>.out` / `.err` to `<outdir>/slurm_jobs/`, and the job runs from there. **The
worker's Unix user must be able to write `--outdir`.**

## When to redeploy

A deployment stores the entrypoint, pool, parameters and job settings at deploy time. Redeploy when
any of these change: a flow moved or renamed (an old entrypoint path makes every run fail), a new
flow or pool, new resources or default parameters, a new branch, or a new venv or ENSCODE path.
Changes *inside* flow code don't need a redeploy: the worker clones the branch for each run.

Schedules (e.g. an hourly `busco-dispatch-automatic`) are set in the Prefect UI, not by the script;
check they are still there after redeploying.
