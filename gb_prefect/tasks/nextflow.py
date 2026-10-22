import dataclasses
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Tuple

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
from gb_prefect.nextflow.metro import MetroRun, register_run
from gb_prefect.nextflow.runner import ignored_errors, run_streaming
from gb_prefect.utils.artifact_utils import (
    artifact_key,
    create_live_map_artifact,
    create_nextflow_run_artifact,
)
from gb_prefect.utils.logging_utils import append_log
from gb_prefect.utils.metro_utils import load_metro_config


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

    Raises NextflowRunError if Nextflow exits non-zero, or exits 0 but ignored failed tasks
    (unless run.fail_on_ignored_errors is False). The artifact shows only this run's part of
    the log, which is appended to across reruns in the same outdir.

    If the run has a metro_map and an nf-metro server is configured (Prefect Variable
    nf-metro-url), the run is registered with it for a live map: -with-weblog is added and a
    link artifact published. Any problem with that is logged as a warning and the run goes on
    without a live map. The server's token is redacted from the command, log and artifacts.
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

    # Prepare output directories and secret strings; note where this run's part of the log starts.
    paths = run_paths(run)
    Path(run.outdir).mkdir(parents=True, exist_ok=True)
    log_start = paths.log.stat().st_size if paths.log.exists() else 0
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

    # Register the run with the nf-metro server for a live map (never stops the run)
    if not dry_run:
        run, secrets = _with_live_map(run, secrets, paths.log, run_options.create_artifact)

    # Build and write the Nextflow command script (as a record: redacted, unlike the argv run)
    command = redact(build_command_script(run, params_file, csv_file), secrets)
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

    # Check this run's output for tasks the pipeline let fail (errorStrategy 'ignore')
    run_log = _read_from(paths.log, log_start)
    ignored = [] if dry_run else ignored_errors(run_log)
    failed = rc != 0 or bool(ignored and run.fail_on_ignored_errors)
    if ignored:
        append_log(
            paths.log,
            f"[{datetime.now()}] {'ERROR' if failed else 'WARNING'}: {len(ignored)} task(s) failed and "
            f"were ignored by the pipeline: {', '.join(ignored)}.\n",
        )
        run_log = _read_from(paths.log, log_start)

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
        ignored_errors=ignored,
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
                "Ignored task failures": len(ignored),
            },
            cmd=command,
            log_text=redact(run_log, secrets),
            rc=rc,
            dry_run=dry_run,
            status=f"FAILED ({len(ignored)} ignored task failure(s))" if failed and rc == 0 else None,
        )

    # Raise an error if the Nextflow run failed
    if failed:
        raise NextflowRunError(result)
    return result


def _with_live_map(
    run: NextflowRun, secrets: List[str], log: Path, create_artifact: bool
) -> Tuple[NextflowRun, List[str]]:
    """Return run with -with-weblog set to its live map's events URL, and secrets with the
    server's token added -- or both unchanged if there is no live map (see _start_live_map)."""
    live_map = _start_live_map(run, log, create_artifact)
    if live_map is None:
        return run, secrets
    metro_run, token = live_map
    return dataclasses.replace(run, weblog_url=metro_run.events_url), secrets + ([token] if token else [])


def _start_live_map(
    run: NextflowRun, log: Path, create_artifact: bool
) -> Optional[Tuple[MetroRun, Optional[str]]]:
    """Register run with the configured nf-metro server; return it and the server's token.

    Returns None -- and the run goes on without a live map -- if the run has no map, already
    has a weblog URL, no server is configured, or anything fails (logged as a warning).
    """
    if not run.metro_map or run.weblog_url:
        return None
    token: Optional[str] = None
    try:
        config = load_metro_config()
        if config is None:
            return None
        token = config.token
        map_text = (Path(run.enscode) / run.metro_map).read_text(encoding="utf-8")
        name = f"{flow_run.name}-{run.name}" if flow_run.name else run.name
        metro_run = register_run(config.url, map_text, name, token)
    except Exception as err:  # pylint: disable=broad-exception-caught  # never stop the run
        _log(
            log,
            f"WARNING: Live map not available ({type(err).__name__}: {err}); running without it.",
            [token],
        )
        return None

    _log(log, f"INFO: Live map: {metro_run.view_url}", [token])
    if create_artifact:
        try:
            create_live_map_artifact(metro_run.view_url, run.name)
        except Exception as err:  # pylint: disable=broad-exception-caught
            _log(log, f"WARNING: Could not publish the live map link ({err}).", [token])
    return metro_run, token


def _log(log: Path, message: str, secrets: List[Optional[str]]) -> None:
    """Append a timestamped message to the log and print it, redacting secrets."""
    message = redact(message, [s for s in secrets if s])
    append_log(log, f"[{datetime.now()}] {message}\n")
    print(message)


def _read_from(path: Path, offset: int) -> str:
    """Return path's content from byte offset on (this run's part of an appended log)."""
    with path.open("rb") as fh:
        fh.seek(offset)
        return fh.read().decode("utf-8", errors="replace")
