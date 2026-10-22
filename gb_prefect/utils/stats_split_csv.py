import csv


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
