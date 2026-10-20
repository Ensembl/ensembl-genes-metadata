from pathlib import Path
from typing import Optional


def run_folder(outdir: str, folder_name: Optional[str], default: str) -> str:
    """Return <outdir>/<folder_name>, or <outdir>/<default> when folder_name is not given.

    folder_name must be a single folder name, not a path.
    """
    if folder_name and "/" in folder_name:
        raise ValueError(f"folder_name must be a single folder name, not a path: {folder_name}")
    return str(Path(outdir) / (folder_name or default))


def check_input_file(path: str, dry_run: bool) -> str:
    """Return path after checking it is an absolute path to an existing file.
    """
    if not Path(path).is_absolute():
        raise ValueError(f"Input file must be an absolute path: {path}")
    if not dry_run and not Path(path).is_file():
        raise FileNotFoundError(f"Input file not found: {path}")
    return path
