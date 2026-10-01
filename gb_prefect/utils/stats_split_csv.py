import csv
from pathlib import Path
from typing import Optional


def write_single_gca_csv(
    gca: str, taxon_id: str, outdir: str, busco_dataset: Optional[str] = None
) -> str:
    """Write a single-row `gca,taxon_id,busco_dataset` CSV for one genome, returning its path."""
    outdir_path = Path(outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)
    csv_path = outdir_path / f"{gca}.csv"

    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["gca", "taxon_id", "busco_dataset"])
        writer.writerow([gca, taxon_id, busco_dataset or ""])

    return str(csv_path)


def read_gca_csv(csv_file: str) -> list[dict[str, str]]:
    """Read a `gca,taxon_id,busco_dataset` CSV into a list of row dicts."""
    with open(csv_file, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        missing = {"gca", "taxon_id"} - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"{csv_file} is missing required column(s): {sorted(missing)}. "
                f"Header found: {reader.fieldnames}"
            )
        return list(reader)


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
