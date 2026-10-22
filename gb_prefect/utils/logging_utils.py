from __future__ import annotations

from pathlib import Path


def append_log(path: Path, msg: str) -> None:
    """Append a message to the log file at path, creating parent directories as needed."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as fh:
        fh.write(msg)
