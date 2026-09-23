import csv
from pathlib import Path


def split_csv(csv_file: str, outdir: str) -> list[str]:
    """Split csv_file into one file per data row (keeping the header), returning the new paths."""
    src = Path(csv_file)
    with open(src, newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = list(reader)

    if len(rows) <= 1:
        return [csv_file]

    split_dir = Path(outdir) / "split_csvs"
    split_dir.mkdir(parents=True, exist_ok=True)

    paths = []
    for i, row in enumerate(rows):
        out_path = split_dir / f"{src.stem}_{i}.csv"
        with open(out_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(header)
            writer.writerow(row)
        paths.append(str(out_path))

    return paths
