import shlex
from datetime import datetime
from pathlib import Path
from typing import Optional

from prefect import task  # type: ignore

from gb_prefect.models.pipeline_options import TaskRunOptions
from gb_prefect.models.script import ScriptResult, ScriptRun, ScriptRunError
from gb_prefect.nextflow.command import redact
from gb_prefect.nextflow.runner import run_streaming
from gb_prefect.utils.artifact_utils import artifact_key, create_nextflow_run_artifact
from gb_prefect.utils.logging_utils import append_log, read_from


@task(log_prints=True)
def run_script(run: ScriptRun, run_options: Optional[TaskRunOptions] = None) -> ScriptResult:
    """Run one command (e.g. a Python script) from run.outdir and return its ScriptResult.

    The script counterpart of run_nextflow_pipeline (tasks/nextflow.py). It runs inside the flow
    run's own Slurm job (slurm-cli work pool), so the script gets that job's resources and the
    pool's setup_commands environment (modules, venv) -- it does not submit a job of its own.

    The command runs directly (no shell), its output streamed line by line to
    <outdir>/script_<name>.log and the Prefect logs. An equivalent bash script is written to
    <outdir>/script_<name>_command.sh as a record you can re-run by hand (secret_env values are
    left out of it). A Prefect artifact shows the command and this run's part of the log.

    On a dry run nothing is executed. Raises ScriptRunError if the command exits non-zero (or
    can't be started: return code 127).
    """
    run_options = run_options or TaskRunOptions()
    outdir = Path(run.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    log = outdir / f"script_{run.name}.log"
    command_file = outdir / f"script_{run.name}_command.sh"
    log_start = log.stat().st_size if log.exists() else 0
    secrets = [value for value in run.secret_env.values() if value]

    command = _command_record(run)
    command_file.write_text(command, encoding="utf-8")
    command_file.chmod(0o755)
    append_log(log, f"[{datetime.now()}] INFO: Command:\n{command}\n")
    print(f"Command for {run.name}:\n{command}")

    if run_options.dry_run:
        rc = 0
        append_log(log, f"[{datetime.now()}] INFO: Dry run enabled; the command was not executed.\n")
    else:
        rc = run_streaming(
            list(run.argv),
            cwd=outdir,
            log_path=log,
            env={**run.env, **run.secret_env},
            secrets=secrets,
            echo=print,
        )
    append_log(log, f"[{datetime.now()}] INFO: Return code {rc}.\n")

    result = ScriptResult(
        name=run.name,
        returncode=rc,
        command=command,
        command_file=str(command_file),
        log_file=str(log),
        outdir=run.outdir,
        dry_run=run_options.dry_run,
    )

    if run_options.create_artifact:
        create_nextflow_run_artifact(
            title=f"Script run {run.name}",
            key=artifact_key(f"script-{run.name}"),
            fields={"Outdir": run.outdir, "Command file": result.command_file, "Log file": result.log_file},
            cmd=command,
            log_text=redact(read_from(log, log_start), secrets),
            rc=rc,
            dry_run=run_options.dry_run,
        )

    if rc != 0:
        raise ScriptRunError(result)
    return result


def _command_record(run: ScriptRun) -> str:
    """The bash equivalent of run, for the record: cd, env exports, then the quoted command.

    secret_env values are not written; their names are listed so a by-hand rerun knows what to set.
    """
    exports = "".join(f"export {name}={shlex.quote(value)}\n" for name, value in run.env.items())
    secret_note = (
        f"# Also set by the task (values not recorded): {', '.join(sorted(run.secret_env))}\n"
        if run.secret_env
        else ""
    )
    return (
        "#!/bin/bash\n"
        "# Record of the command run by gb_prefect's run_script (executed directly, not via this file).\n"
        "set -euo pipefail\n"
        f"cd {shlex.quote(run.outdir)}\n"
        f"{exports}{secret_note}"
        f"{shlex.join(run.argv)}\n"
    )
