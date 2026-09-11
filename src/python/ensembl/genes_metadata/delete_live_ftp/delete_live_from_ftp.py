"""Safely delete GCA directories listed in a reviewed manifest."""
import argparse
import re
import shutil
from pathlib import Path

PRE_RELEASE_ROOT = Path("/nfs/ftp/public/databases/ensembl/pre-release").resolve()
GCA_DIRECTORY = re.compile(r"^GCA_\d+\.\d+$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Delete GCA directories from the pre-release FTP tree."
    )
    parser.add_argument("manifest", type=Path)
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Actually delete directories; otherwise only validate and report.",
    )
    parser.add_argument("--root", type=Path, default=PRE_RELEASE_ROOT)
    return parser.parse_args()


def read_manifest(manifest: Path) -> list[Path]:
    paths = []
    with manifest.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("path\t"):
                continue
            fields = line.split("\t")
            if len(fields) >= 5:
                if fields[4] == "delete":
                    paths.append(Path(fields[0]))
            elif len(fields) == 1:
                paths.append(Path(fields[0]))
            else:
                raise ValueError(f"Unrecognised manifest row: {line}")
    return paths


def validate_path(path: Path, root: Path) -> Path:
    resolved = path.resolve(strict=False)
    root = root.resolve()
    try:
        relative = resolved.relative_to(root)
    except ValueError as exc:
        raise ValueError(f"Path is outside pre-release root: {path}") from exc

    if len(relative.parts) != 2 or not GCA_DIRECTORY.fullmatch(relative.parts[1]):
        raise ValueError(f"Path is not a species/GCA directory: {path}")
    if path.is_symlink():
        raise ValueError(f"Refusing to delete symlink: {path}")
    return resolved


def main() -> None:
    args = parse_args()
    if not args.manifest.exists():
        raise SystemExit(f"File does not exist: {args.manifest}")

    paths = read_manifest(args.manifest)
    validated = [validate_path(path, args.root) for path in paths]
    mode = "EXECUTE" if args.execute else "DRY RUN"
    print(f"{mode}: validated {len(validated)} manifest paths")

    for gca_path in validated:
        if not gca_path.exists() or not gca_path.is_dir():
            print(f"SKIP: path does not exist or is not a directory: {gca_path}")
            continue
        if not args.execute:
            print(f"WOULD DELETE: {gca_path}")
            continue

        print(f"DELETING: {gca_path}")
        shutil.rmtree(gca_path)
        parent = gca_path.parent
        if parent.exists() and not any(parent.iterdir()):
            print(f"DELETING EMPTY SPECIES DIRECTORY: {parent}")
            parent.rmdir()


if __name__ == "__main__":
    main()
