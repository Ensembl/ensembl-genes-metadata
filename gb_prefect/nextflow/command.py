"""Turn a NextflowRun into what launching it needs -- plain Python, no Prefect, nothing executed here.

- nextflow_argv / nextflow_env: the argument list and extra environment the runner executes
  (directly, no shell).
- build_command_script: the same command as a bash script, written next to the run as a record.
- write_params_file / write_csv_input: the owner-only -params-file and the generated input CSV.
- secret_strings / redact: what to hide from logs and artifacts.
"""

import csv
import json
import os
import re
import shlex
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional

from gb_prefect.models.nextflow import CsvInput, NextflowRun

REDACTED = "<REDACTED>"

# Keys whose values are treated as secrets wherever they appear inside secret_params: the whole
# key or its last "_"-separated part, e.g. password, pass, db_password, slack_token, api_key.
_SECRET_KEY = re.compile(r"(^|_)(pass|passwd|password|pwd|secret|token|api_?key)$", re.IGNORECASE)

# Shorter strings are never redacted: replacing them everywhere would garble the logs.
MIN_SECRET_LENGTH = 4


@dataclass(frozen=True)
class RunPaths:
    """Files of one run, all inside its outdir and named after it."""

    log: Path
    command_file: Path
    params_file: Path
    csv_file: Path


def run_paths(run: NextflowRun) -> RunPaths:
    """Return the paths of run's log, command script, params file and input CSV."""
    outdir = Path(run.outdir)
    return RunPaths(
        log=outdir / f"nextflow_{run.name}.log",
        command_file=outdir / f"nextflow_{run.name}_command.sh",
        params_file=outdir / f"nextflow_{run.name}_params.json",
        csv_file=outdir / f"{run.name}.csv",
    )


def _render_param(key: str, value: Any) -> List[str]:
    """Render one pipeline parameter as command-line words (empty if it's left out)."""
    if value is None:
        return []
    if value is True:
        return [f"--{key}"]
    if value is False:
        return [f"--{key}", "false"]
    if isinstance(value, (str, int, float, Path)):
        return [f"--{key}", str(value)]
    raise TypeError(
        f"Parameter {key!r} has unsupported type {type(value).__name__}; "
        "pass structured values through secret_params (a -params-file) instead."
    )


def nextflow_args(run: NextflowRun, params_file: Optional[Path], csv_file: Optional[Path]) -> List[List[str]]:
    """Return the nextflow command as groups of words, one group per rendered line.

    params_file and csv_file are passed in rather than derived, so the caller decides whether
    they exist (e.g. no params file when secret_params is empty).
    """
    head = ["nextflow"]
    if run.config:
        head += ["-C", f"{run.enscode}/{run.config}"]
    head += ["run", f"{run.enscode}/{run.pipeline}"]

    groups: List[List[str]] = [head]
    if csv_file is not None:
        if run.csv_input is None:
            raise ValueError("csv_file given but run has no csv_input")
        groups.append(_render_param(run.csv_input.param, str(csv_file)))
    for key, value in run.params.items():
        words = _render_param(key, value)
        if words:
            groups.append(words)
    if params_file is not None:
        groups.append(["-params-file", str(params_file)])
    if run.profile:
        groups.append(["-profile", run.profile])
    if run.with_report:
        groups.append(["-with-report"])
    if run.dag_file:
        groups.append(["-with-dag", str(Path(run.outdir) / run.dag_file)])
    if run.weblog_url:
        groups.append(["-with-weblog", run.weblog_url])
    if run.resume:
        groups.append(["-resume"])
    # Plain line-by-line progress instead of the interactive display: the output is a log file.
    groups.append(["-ansi-log", "false"])
    return groups


def nextflow_argv(run: NextflowRun, params_file: Optional[Path], csv_file: Optional[Path]) -> List[str]:
    """Return the nextflow command as the flat argument list that is executed."""
    return [word for group in nextflow_args(run, params_file, csv_file) for word in group]


def nextflow_env(run: NextflowRun) -> Dict[str, str]:
    """Return the environment variables a run sets on top of the inherited environment."""
    return {"ENSCODE": run.enscode, "NXF_DISABLE_CHECK_LATEST": "true"}


def build_command_script(run: NextflowRun, params_file: Optional[Path], csv_file: Optional[Path]) -> str:
    """Return a bash script equivalent to the run, kept as a record next to its outputs.

    The task does not execute this script -- it runs nextflow_argv directly from run.outdir with
    nextflow_env. The script renders the same directory, environment and arguments, so it can
    be read in the log and artifact, or re-run by hand (`bash <command file>`) for debugging.

    Every value is shell-quoted. The script never contains secret_params -- only the path of
    the params file holding them.
    """
    groups = nextflow_args(run, params_file, csv_file)
    command = " \\\n    ".join(" ".join(shlex.quote(word) for word in group) for group in groups)
    exports = "\n".join(f"export {name}={shlex.quote(value)}" for name, value in nextflow_env(run).items())
    return f"""#!/bin/bash
# Record of the Nextflow run launched by gb_prefect (executed directly, not via this file).
set -euo pipefail
cd {shlex.quote(run.outdir)}
{exports}
{command}
"""


def write_params_file(path: Path, secret_params: Mapping[str, Any]) -> None:
    """Write secret_params as JSON to path, readable by the owner only.

    The file is created with mode 0600, so its content is never readable by anyone else,
    even briefly; an existing file is truncated and set back to 0600.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
    os.fchmod(fd, 0o600)
    with os.fdopen(fd, "w", encoding="utf-8") as fh:
        json.dump(dict(secret_params), fh)


def write_csv_input(path: Path, csv_input: CsvInput) -> None:
    """Write csv_input's header and rows to path."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(csv_input.columns), restval="", extrasaction="raise")
        writer.writeheader()
        for row in csv_input.rows:
            writer.writerow({key: "" if value is None else value for key, value in row.items()})


def secret_strings(secret_params: Mapping[str, Any]) -> List[str]:
    """Return the strings to redact from logs for secret_params.

    - every top-level string value (e.g. a whole metadata_params_string);
    - every top-level dict or list, as JSON (e.g. a dumped asm_metadata);
    - any nested value under a secret-looking key (see _SECRET_KEY), including inside
      top-level strings that hold a JSON object.

    Top-level numbers and non-secret nested values (host, port, user) are not collected, and
    nothing shorter than MIN_SECRET_LENGTH is -- either would replace common text all over the
    logs. Longest first, so a JSON blob is replaced whole before a password inside it.
    """
    found = set()

    def visit(value: Any) -> None:
        if isinstance(value, Mapping):
            for key, inner in value.items():
                is_scalar = isinstance(inner, (str, int)) and not isinstance(inner, bool)
                if is_scalar and _SECRET_KEY.search(str(key)):
                    found.add(str(inner))
                visit(inner)
        elif isinstance(value, (list, tuple)):
            for inner in value:
                visit(inner)

    for value in secret_params.values():
        if isinstance(value, str):
            found.add(value)
            try:
                visit(json.loads(value))
            except ValueError:
                pass
        elif isinstance(value, (Mapping, list, tuple)):
            found.add(json.dumps(value))
            visit(value)
    return sorted((s for s in found if len(s) >= MIN_SECRET_LENGTH), key=len, reverse=True)


def redact(text: str, secrets: List[str]) -> str:
    """Replace every string in secrets found in text with REDACTED."""
    for secret in secrets:
        text = text.replace(secret, REDACTED)
    return text
