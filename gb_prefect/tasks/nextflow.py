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
    """Universal Nextflow pipeline runner that takes a NextflowRun object, executes it, 
    and returns a NextflowResult

    run: NextflowRun object containing the pipeline run configuration.
    run_options: Optional TaskRunOptions object containing additional run options.
    resume_on_retry: Boolean flag indicating whether to resume the run on retry.
    """
    # Initialize run options and determine if this is a dry run.
    run_options = run_options or TaskRunOptions()
    dry_run = run_options.dry_run

    # Determine if the run should be resumed based on the retry count and the resume_on_retry flag.
    if resume_on_retry and not run.resume and (flow_run.run_count or 0) > 1:
        run = dataclasses.replace(run, resume=True)

    # Generate DAG file if specified
    if run.dag_file:
        run = dataclasses.replace(run, dag_file=datetime.now().strftime(run.dag_file))

    # Prepare output directories and secret strings.
    paths = run_paths(run)
    Path(run.outdir).mkdir(parents=True, exist_ok=True)
    secrets = secret_strings(run.secret_params)
    files = {}

    # Initialize CSV input file
    csv_file = None
    if run.csv_input is not None:
        write_csv_input(paths.csv_file, run.csv_input)
        csv_file = paths.csv_file
        files["csv"] = str(csv_file)
        append_log(paths.log, f"[{datetime.now()}] INFO: Wrote input CSV {csv_file}.\n")

    # Initialize params file
    params_file = None
    if run.secret_params:
        params_file = paths.params_file
        files["params"] = str(params_file)
        if not dry_run:
            write_params_file(params_file, run.secret_params)
            append_log(paths.log, f"[{datetime.now()}] INFO: Wrote params file {params_file}.\n")

    # Build and write the Nextflow command script
    command = build_command_script(run, params_file, csv_file)
    paths.command_file.write_text(command, encoding="utf-8")
    paths.command_file.chmod(0o755)
    append_log(paths.log, f"[{datetime.now()}] INFO: Nextflow command:\n{command}\n")
    print(f"Nextflow command for {run.name}:\n{command}")

    # Execute the Nextflow command or perform a dry run
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

    # collect the results into a NextflowResult object
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

    # Create an artifact for the Nextflow run if requested
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
    
    # Raise an error if the Nextflow run failed
    if rc != 0:
        raise NextflowRunError(result)
    return result
