from __future__ import annotations

import re
from pathlib import Path

IGNORE = re.compile(r"Cannot infer \$ENSEMBL_VERSION from module's version")


def append_log(path: Path, msg: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as fh:
        fh.write(msg)


def filter_stderr(s: str) -> str:
    if not s:
        return ""
    return "\n".join(line for line in s.splitlines() if not IGNORE.search(line)).strip()
