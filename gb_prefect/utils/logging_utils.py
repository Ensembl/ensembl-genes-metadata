from __future__ import annotations

import re
from pathlib import Path

IGNORE = re.compile(r"Cannot infer \$ENSEMBL_VERSION from module's version")


def append_log(path: Path, msg: str) -> None:
    """Append a message to the log file at path, creating parent directories as needed."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as fh:
        fh.write(msg)


def filter_stderr(s: str) -> str:
    """Strip out known-noisy lines (e.g. module version warnings) from stderr output."""
    if not s:
        return ""
    return "\n".join(line for line in s.splitlines() if not IGNORE.search(line)).strip()
