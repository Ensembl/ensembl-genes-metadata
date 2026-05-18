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
