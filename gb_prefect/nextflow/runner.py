"""Run a command, streaming its output to a log file as it goes."""

import os
import re
import shlex
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Callable, List, Mapping, Optional

from gb_prefect.nextflow.command import redact

# Return code when the command can't be started at all (as a shell reports "command not found").
NOT_STARTED_RC = 127

# How long a stopped command gets to shut down after SIGTERM before it is killed. Nextflow uses
# this time to cancel the Slurm jobs it submitted.
STOP_TIMEOUT_S = 60

# Terminal colour/formatting escape codes, e.g. in nf-schema's parameter summary, which is coloured
# even with -ansi-log false.
_ANSI_ESCAPE = re.compile(r"\x1b\[[0-?]*[ -/]*[@-~]")

# The line Nextflow prints when a task fails for good and the pipeline's errorStrategy is 'ignore':
# NOTE: Process `RUN_BUSCO:FETCH_GENOME (GCA_1.1:genome)` terminated with an error exit status (1)
# -- Error is ignored. A failure that is retried prints "-- Execution is retried" instead.
_IGNORED_ERROR = re.compile(
    r"Process `([^`]+)` terminated with an error exit status \((-?\d+)\) -- Error is ignored"
)


def strip_ansi(text: str) -> str:
    """Remove terminal colour/formatting escape codes from text."""
    return _ANSI_ESCAPE.sub("", text)


def ignored_errors(output: str) -> List[str]:
    """Return the tasks Nextflow's output says failed and were ignored, as "<task> (exit <status>)"."""
    return [f"{task} (exit {status})" for task, status in _IGNORED_ERROR.findall(output)]


def run_streaming(  # pylint: disable=too-many-arguments
    argv: List[str],
    cwd: Path,
    log_path: Path,
    env: Optional[Mapping[str, str]] = None,
    secrets: Optional[List[str]] = None,
    echo: Optional[Callable[[str], None]] = None,
) -> int:
    """Run argv from cwd -- no shell -- and return its exit code.

    env is added on top of the inherited environment: nextflow and the venv come from the
    process running this, i.e. the slurm-cli worker's setup_commands.

    stdout and stderr are merged and each line is appended to log_path as it arrives
    (earlier log content is kept), and passed to echo -- e.g. print, which a Prefect task
    with log_prints=True turns into live run logs. Lines have terminal escape codes removed and
    are redacted with secrets first.
    Output is read as UTF-8; undecodable bytes are replaced rather than ending the run.

    If argv can't be started (e.g. nextflow not on PATH), that is logged and NOT_STARTED_RC
    returned, so the caller fails the run through its normal path. If anything goes wrong
    while reading (an exception, or the caller being cancelled), the command is stopped --
    SIGTERM, then SIGKILL after STOP_TIMEOUT_S -- before the error is re-raised, so it is
    never left running unattended.
    """
    secrets = secrets or []
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as log:

        def emit(line: str) -> None:
            log.write(line if line.endswith("\n") else f"{line}\n")
            log.flush()
            if echo:
                echo(line.rstrip("\n"))

        log.write(f"[{datetime.now()}] INFO: Running in {cwd}: {redact(shlex.join(argv), secrets)}\n")
        log.flush()
        try:
            proc = subprocess.Popen(  # pylint: disable=consider-using-with  # closed in finally
                argv,
                cwd=cwd,
                env={**os.environ, **(env or {})},
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                encoding="utf-8",
                errors="replace",
                bufsize=1,
            )
        except OSError as err:
            emit(
                f"[{datetime.now()}] ERROR: Could not start {argv[0]}: {err}. "
                "If it is not on PATH, check the work pool's setup_commands (module load, venv).\n"
            )
            return NOT_STARTED_RC

        try:
            assert proc.stdout is not None
            for line in proc.stdout:
                emit(redact(strip_ansi(line), secrets))
            return proc.wait()
        except BaseException:
            _stop(proc)
            raise
        finally:
            if proc.stdout is not None:
                proc.stdout.close()


def _stop(proc: subprocess.Popen) -> None:
    """Terminate proc, giving it STOP_TIMEOUT_S to exit cleanly before killing it."""
    if proc.poll() is not None:
        return
    proc.terminate()
    try:
        proc.wait(timeout=STOP_TIMEOUT_S)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.wait()
