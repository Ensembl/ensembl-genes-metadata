from __future__ import annotations

import subprocess
from pathlib import Path


def run_cmd_bash_capture(
    cmd: str,
    log_path: Path | None = None,
) -> subprocess.CompletedProcess[str]:
    """Run cmd in a login bash shell, capturing stdout/stderr and optionally logging them."""
    full = f"source /etc/bashrc && {cmd}"
    proc = subprocess.run(
        ["bash", "-lc", full],
        text=True,
        capture_output=True,
        check=False,
    )

    if log_path:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text((proc.stdout or "") + "\n" + (proc.stderr or ""), encoding="utf-8")

    return proc
