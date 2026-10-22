import re
from typing import Any, Dict, Optional

from prefect.artifacts import create_markdown_artifact  # type: ignore


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
    status: Optional[str] = None,
) -> None:
    """Publish a Prefect markdown artifact summarizing any Nextflow run.

    fields are shown as a table under the status; key must be lowercase letters, digits and
    dashes (see artifact_key). status overrides the status derived from rc.
    """
    status = status or ("SUCCESS" if rc == 0 else f"FAILED (rc={rc})")
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
