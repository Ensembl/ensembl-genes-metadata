from typing import Optional
import os


def resolve_enscode(enscode: Optional[str], dry_run: bool) -> str:
    """Resolve ENSCODE from the given value, the environment, or a dry-run placeholder."""
    resolved = enscode or os.environ.get("ENSCODE")
    if not resolved:
        if dry_run:
            return "<ENSCODE>"
        raise ValueError(
            "ENSCODE is required when dry_run=False. Pass enscode=..., set ENSCODE, or run with dry_run=True."
        )
    return resolved
