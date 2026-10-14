import re
from typing import Any, Dict

from prefect.artifacts import create_markdown_artifact  # type: ignore


def create_busco_run_artifact(
    csv_file: str,
    outdir: str,
    command_file: str,
    cmd: str,
    log_text: str,
    rc: int,
    dry_run: bool,
) -> None:
    """Publish a Prefect markdown artifact summarizing a BUSCO Nextflow run."""
    status = "SUCCESS" if rc == 0 else f"FAILED (rc={rc})"
    dry_run_badge = " *(dry run)*" if dry_run else ""
    markdown = f"""\
# BUSCO run log

**Status:** {status}{dry_run_badge}

| Field | Value |
|---|---|
| CSV file | `{csv_file}` |
| Outdir | `{outdir}` |
| Command file | `{command_file}` |
| Return code | `{rc}` |

## Nextflow command

```bash
{cmd.strip()}
```

## Log output

```
{log_text.strip()}
```
"""
    create_markdown_artifact(
        key="busco-run-log",
        description="Nextflow BUSCO run log",
        markdown=markdown,
    )


def create_registry_run_artifact(
    date: str,
    outdir: str,
    command_file: str,
    cmd: str,
    log_text: str,
    rc: int,
    dry_run: bool,
) -> None:
    """Publish a Prefect markdown artifact summarizing an assembly registry Nextflow run."""
    status = "SUCCESS" if rc == 0 else f"FAILED (rc={rc})"
    dry_run_badge = " *(dry run)*" if dry_run else ""
    markdown = f"""\
# Assembly registry run log

**Status:** {status}{dry_run_badge}

| Field | Value |
|---|---|
| Date | `{date}` |
| Outdir | `{outdir}` |
| Command file | `{command_file}` |
| Return code | `{rc}` |

## Nextflow command

```bash
{cmd.strip()}
```

## Log output

```
{log_text.strip()}
```
"""
    create_markdown_artifact(
        key="registry-run-log",
        description="Nextflow Assembly registry run log",
        markdown=markdown,
    )


# The artifact shows the end of the log; the full log stays in the run's outdir.
ARTIFACT_LOG_TAIL_LINES = 300


def create_nextflow_run_artifact(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    title: str,
    key: str,
    fields: Dict[str, Any],
    cmd: str,
    log_text: str,
    rc: int,
    dry_run: bool,
) -> None:
    """Publish a Prefect markdown artifact summarizing any Nextflow run.

    fields are shown as a table under the status; key must be lowercase letters, digits and
    dashes (see artifact_key).
    """
    status = "SUCCESS" if rc == 0 else f"FAILED (rc={rc})"
    dry_run_badge = " *(dry run)*" if dry_run else ""
    rows = "\n".join(f"| {name} | `{value}` |" for name, value in {**fields, "Return code": rc}.items())
    log_lines = log_text.strip().splitlines()
    if len(log_lines) > ARTIFACT_LOG_TAIL_LINES:
        hidden = len(log_lines) - ARTIFACT_LOG_TAIL_LINES
        log_lines = [f"... ({hidden} earlier lines in the log file)", *log_lines[-ARTIFACT_LOG_TAIL_LINES:]]
    log_block = "\n".join(log_lines)
    markdown = f"""\
# {title}

**Status:** {status}{dry_run_badge}

| Field | Value |
|---|---|
{rows}

## Nextflow command

```bash
{cmd.strip()}
```

## Log output

```
{log_block}
```
"""
    create_markdown_artifact(key=key, description=title, markdown=markdown)


def artifact_key(name: str) -> str:
    """Turn name into a valid Prefect artifact key (lowercase letters, digits and dashes)."""
    return re.sub(r"-+", "-", re.sub(r"[^a-z0-9-]", "-", name.lower())).strip("-") or "nextflow-run"
