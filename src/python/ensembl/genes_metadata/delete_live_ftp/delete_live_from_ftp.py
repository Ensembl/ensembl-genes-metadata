import shutil
from pathlib import Path
import sys

if len(sys.argv) != 2:
    print(f"Usage: {sys.argv[0]} <paths_file>")
    sys.exit(1)

paths_file = Path(sys.argv[1])
if not paths_file.exists():
    print(f"File does not exist: {paths_file}")
    sys.exit(1)

with open(paths_file) as f:
    paths = [line.strip() for line in f if line.strip()]

for path_str in paths:
    gca_path = Path(path_str)
    if gca_path.exists() and gca_path.is_dir():
        print(f"Deleting GCA folder: {gca_path}")
        shutil.rmtree(gca_path)

        parent = gca_path.parent
        if parent.exists() and not any(parent.iterdir()):
            print(f"Parent folder is empty. Deleting parent: {parent}")
            parent.rmdir()
    else:
        print(f"Path does not exist or is not a directory: {gca_path}")
