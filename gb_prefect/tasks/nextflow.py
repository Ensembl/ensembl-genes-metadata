import dataclasses
from datetime import datetime
from pathlib import Path
from typing import Optional

from prefect import task  # type: ignore
from prefect.runtime import flow_run  # type: ignore

from gb_prefect.models.nextflow import NextflowResult, NextflowRun, NextflowRunError
from gb_prefect.models.pipeline_options import TaskRunOptions
from gb_prefect.nextflow.command import (
    build_command_script,
    nextflow_argv,
    nextflow_env,
    redact,
    run_paths,
    secret_strings,
    write_csv_input,
    write_params_file,
)
from gb_prefect.nextflow.runner import run_streaming
from gb_prefect.utils.artifact_utils import artifact_key, create_nextflow_run_artifact
from gb_prefect.utils.logging_utils import append_log


@task(log_prints=True)
def run_nextflow_pipeline(
    run: NextflowRun,
    run_options: Optional[TaskRunOptions] = None,
    resume_on_retry: bool = True,
) -> NextflowResult:
    """Write run's files to its outdir and run Nextflow there; the one task every Nextflow flow uses.

    Assumes it executes inside the Slurm job the slurm-cli worker submitted for the flow run,
    with nextflow on PATH via the worker's setup_commands -- it runs `nextflow run` directly.

    Steps: write the input CSV (if any) and the owner-only -params-file (if any secret_params),
    run nextflow from outdir -- directly, no shell -- streaming its output to the log and the
    Prefect logs, then publish an artifact (run_options.create_artifact). An equivalent bash
    script is written to the command file as a record, re-runnable by hand; it isn't executed.

    resume_on_retry adds -resume when the flow run is a retry (run_count > 1): its outdir is
    the same, so Nextflow picks up the cached tasks.

    On a dry run nothing is executed and the params file is not written (the command still
    refers to it). A non-zero exit raises NextflowRunError, so every Nextflow flow fails the
    same way -- returning a Failed state carrying a dict crashes Prefect's result resolution.
    """
    run_options = run_options or TaskRunOptions()
    dry_run = run_options.dry_run
    if resume_on_retry and not run.resume and (flow_run.run_count or 0) > 1:
        run = dataclasses.replace(run, resume=True)

    paths = run_paths(run)
    Path(run.outdir).mkdir(parents=True, exist_ok=True)
    secrets = secret_strings(run.secret_params)
    files = {}

    csv_file = None
    if run.csv_input is not None:
        write_csv_input(paths.csv_file, run.csv_input)
        csv_file = paths.csv_file
        files["csv"] = str(csv_file)
        append_log(paths.log, f"[{datetime.now()}] INFO: Wrote input CSV {csv_file}.\n")

    params_file = None
    if run.secret_params:
        params_file = paths.params_file
        files["params"] = str(params_file)
        if not dry_run:
            write_params_file(params_file, run.secret_params)
            append_log(paths.log, f"[{datetime.now()}] INFO: Wrote params file {params_file}.\n")

    command = build_command_script(run, params_file, csv_file)
    paths.command_file.write_text(command, encoding="utf-8")
    paths.command_file.chmod(0o755)
    append_log(paths.log, f"[{datetime.now()}] INFO: Nextflow command:\n{command}\n")
    print(f"Nextflow command for {run.name}:\n{command}")

    if dry_run:
        rc = 0
        append_log(paths.log, f"[{datetime.now()}] INFO: Dry run enabled; Nextflow was not executed.\n")
    else:
        rc = run_streaming(
            nextflow_argv(run, params_file, csv_file),
            cwd=Path(run.outdir),
            log_path=paths.log,
            env=nextflow_env(run),
            secrets=secrets,
            echo=print,
        )
    append_log(paths.log, f"[{datetime.now()}] INFO: Return code {rc}.\n")

    result = NextflowResult(
        name=run.name,
        returncode=rc,
        command=command,
        command_file=str(paths.command_file),
        log_file=str(paths.log),
        outdir=run.outdir,
        dry_run=dry_run,
        resumed=run.resume,
        files=files,
    )

    if run_options.create_artifact:
        create_nextflow_run_artifact(
            title=f"Nextflow run {run.name}",
            key=artifact_key(f"nextflow-{run.name}"),
            fields={
                "Pipeline": run.pipeline,
                "Outdir": run.outdir,
                "Command file": result.command_file,
                "Log file": result.log_file,
                "Resumed": run.resume,
            },
            cmd=command,
            log_text=redact(paths.log.read_text(encoding="utf-8"), secrets),
            rc=rc,
            dry_run=dry_run,
        )

    if rc != 0:
        raise NextflowRunError(result)
    return result
