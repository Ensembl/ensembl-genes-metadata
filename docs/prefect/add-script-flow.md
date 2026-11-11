# Add a script flow

There are two ways to run Python from a flow. Pick by what the work is:

| | **Script flow** (`run_script`) | **Python flow** (in-process) |
|---|---|---|
| Runs | a command: a Python script, a binary, anything with a command line | Python functions, imported and called by the flow |
| Where | inside the flow run's own Slurm job (`slurm-cli` pool) | in the flow run's process (usually the `genebuild-dispatch` process pool) |
| Good for | real work: long, heavy, or a script that already has a CLI | quick work: DB queries, API calls, triggering other deployments |
| You get | live logs, a re-runnable command record, a run-summary artifact, `ScriptRunError` on failure | Prefect's task logging; you write the steps |
| Example | [`gb_is_reference.py`](https://github.com/Ensembl/ensembl-genes-metadata/blob/main/gb_prefect/flows/gb_is_reference.py) | [`gb_busco_dispatch_automatic.py`](https://github.com/Ensembl/ensembl-genes-metadata/blob/main/gb_prefect/flows/statistics/gb_busco_dispatch_automatic.py) |

Read the [overview](index.md) first if Prefect is new to you.

## Script flow with `run_script`

The example adds a flow for a made-up script, `src/python/species_report.py`, which takes a GCA list
and writes a report. The flow builds a
[`ScriptRun`](https://github.com/Ensembl/ensembl-genes-metadata/blob/main/gb_prefect/models/script.py)
and passes it to
[`run_script`](https://github.com/Ensembl/ensembl-genes-metadata/blob/main/gb_prefect/tasks/script.py).

### 1. Write the flow

`gb_prefect/flows/assembly_metadata/gb_species_report.py`:

```python
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional

from prefect import flow  # type: ignore

from gb_prefect.models.pipeline_options import TaskRunOptions
from gb_prefect.models.script import ScriptRun
from gb_prefect.tasks.script import run_script
from gb_prefect.utils.credentials_utils import PipelineCredentials, resolve_credentials
from gb_prefect.utils.enscode_utils import resolve_enscode
from gb_prefect.utils.path_utils import check_input_file, run_folder


@flow(name="gb_species_report", log_prints=True)
def gb_species_report_flow(
    gca_list: str,
    outdir: str,
    folder_name: Optional[str] = None,
    enscode: Optional[str] = None,
    dry_run: bool = False,
    credentials: Optional[PipelineCredentials] = None,
    metadata_secret_block: Optional[str] = None,
) -> str:
    """Write a species report for the GCA accessions in gca_list; return the report's path.

    gca_list: absolute path to a file with one GCA accession per line.
    outdir, folder_name: the run goes to <outdir>/<folder_name> (default species_report_<today>).
    enscode: ENSCODE directory holding ensembl-genes-metadata; defaults to $ENSCODE.
    dry_run: write and log the command without running it or loading credentials.
    credentials, metadata_secret_block: DB credentials, explicit or from the Secret block.
    """
    run_outdir = run_folder(outdir, folder_name, f"species_report_{datetime.now().strftime('%Y-%m-%d')}")
    report = Path(run_outdir) / "species_report.tsv"
    credentials = resolve_credentials(credentials, dry_run=dry_run, metadata_secret_block=metadata_secret_block)

    run = ScriptRun(
        name="species_report",
        argv=[
            sys.executable,                     # the venv the work pool activates
            f"{resolve_enscode(enscode, dry_run)}/ensembl-genes-metadata/src/python/species_report.py",
            "--gca-list", check_input_file(gca_list, dry_run),
            "--output", str(report),
        ],
        outdir=run_outdir,
        secret_env={"METADATA_PARAMS": credentials.metadata_params_string},
    )
    run_script(run, TaskRunOptions(dry_run=dry_run))
    return str(report)
```

`ScriptRun` fields:

| Field | Notes |
|---|---|
| `name` | Names `script_<name>.log`, `script_<name>_command.sh` and the artifact; unique within `outdir`. |
| `argv` | The command as a list. It runs **without a shell**: no quoting needed, and `$VARS`, `*`, pipes and `;` are passed literally, not interpreted. For a pipeline of commands, write a script. |
| `outdir` | The command runs from here; the log and record are written here. |
| `env` | Extra environment variables (on top of the pool's `setup_commands` environment). |
| `secret_env` | Environment variables holding secrets. Their values are redacted from the log and artifact and never written to the command record. **Credentials go here, never in `argv`**: command lines show up in process lists and Nextflow-style logs. The script reads them from its environment, e.g. `json.loads(os.environ["METADATA_PARAMS"])`. |

`sys.executable` is the Python of the venv the work pool's `setup_commands` activate, so the script
runs with the same packages as the flow. The script comes from **ENSCODE**, like pipelines do: the
version that runs is the ENSCODE checkout's.

`run_script` raises `ScriptRunError` if the command exits non-zero (127 if it can't be started), which
fails the flow run. On a dry run nothing is executed; the command record still shows what would run.

### 2. Add the deployment

In
[`deploy_flows.py`](https://github.com/Ensembl/ensembl-genes-metadata/blob/main/gb_prefect/setup/deploy_flows.py)'s
`FLOWS`, on a `slurm-cli` pool, with the resources the **script** needs (it runs in the flow's job):

```python
    "species_report": FlowDeployment(
        name="gb_species_report",
        entrypoint="gb_prefect/flows/assembly_metadata/gb_species_report.py:gb_species_report_flow",
        work_pool="metadata",
        description="Species report",
        tags=["metadata_report"],
        memory=2,
        time_limit=1,
    ),
```

Leave `runs_nextflow` at its default (`True`): despite the name it means "a flow that runs as a Slurm
job, with `outdir` and `enscode` parameters".

### 3. Test, deploy and run

As for a Nextflow flow ([Testing](testing.md), step 6 of [Add a Nextflow flow](add-nextflow-flow.md)):
a dry run first, then a `_test` deployment. The run folder holds `script_species_report.log` (the
command and its output) and `script_species_report_command.sh` (re-run it by hand with `bash`, after
setting the variables it lists as not recorded).

## Python flow (in-process)

When the work is a few Python calls, skip the subprocess: import the functions and call them from
tasks. This is how the BUSCO dispatchers work.

```python
from typing import Optional

from prefect import flow, task  # type: ignore

from gb_metadata.update_busco_events import get_pending_candidates  # type: ignore
from gb_prefect.utils.credentials_utils import PipelineCredentials, resolve_credentials


@task(log_prints=True)
def count_pending(db_params: dict, limit: int) -> int:
    """Return how many GCAs are waiting for BUSCO (up to limit)."""
    pending = get_pending_candidates(db_params, limit)
    print(f"{len(pending)} GCA(s) pending")
    return len(pending)


@flow(name="gb_busco_pending", log_prints=True)
def gb_busco_pending_flow(
    limit: int = 1000,
    credentials: Optional[PipelineCredentials] = None,
    metadata_secret_block: Optional[str] = None,
) -> int:
    """Report how many GCAs are waiting for BUSCO."""
    credentials = resolve_credentials(credentials, metadata_secret_block=metadata_secret_block)
    return count_pending(credentials.metadata_db_params, limit)
```

- Wrap each meaningful step in a `@task`: it shows as a step in the UI, with its own logs and state.
  `print` goes to the Prefect logs because of `log_prints=True`.
- `credentials.metadata_db_params` is the DB connection as a dict (for `pymysql`);
  `metadata_params_string` is the same as JSON, for pipelines and scripts.
- Deploy it on the process pool, with no Slurm settings and only its own parameters:

```python
    "busco_pending": FlowDeployment(
        name="gb-busco-pending",
        entrypoint="gb_prefect/flows/statistics/gb_busco_pending.py:gb_busco_pending_flow",
        work_pool="dispatch",
        runs_nextflow=False,
        description="BUSCO pending count",
        tags=["busco"],
    ),
```

!!! warning "Where in-process code comes from"
    The flow's own code (`gb_prefect`) is cloned from the deployment's branch for each run, but
    `gb_metadata` (`src/python/`) is imported from the **installed package in the venv**. After
    changing `src/python/`, reinstall it in the venv (`pip install -e <checkout>` and `git pull`
    there), or the flow runs the old code or fails with an `ImportError`. Script flows don't have
    this problem: they run the script from ENSCODE.

Keep in-process flows quick: they hold a process on the worker's machine for as long as they run.
Anything long or heavy belongs in a script flow on a `slurm-cli` pool.
