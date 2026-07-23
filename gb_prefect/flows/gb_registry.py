import argparse
import sys
from datetime import datetime
from pathlib import Path

from prefect import flow  # type: ignore
from gb_prefect.tasks.registry import register_assemblies  # pylint: disable=wrong-import-position


@flow(name="gb_registry", log_prints=True)
def gb_registry_flow(
    date: str,
    outdir: str,
    enscode: str,
    asm_venv: str,
    dry_run: bool = False,
):
    """Run the assembly registry Nextflow pipeline for the given date."""
    return register_assemblies(
        date=date,
        outdir=f"{outdir}/{date}",
        enscode=enscode,
        asm_venv=asm_venv,
        dry_run=dry_run,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dry-run", action="store_true", help="Create the Nextflow command without running it."
    )
    parser.add_argument("--date", required=True, help="Date for the registry run (e.g., MM-DD-YYYY).")
    parser.add_argument("--outdir", required=True, help="Base output directory.")
    parser.add_argument("--enscode", required=True, help="Path to ENSCODE directory.")
    parser.add_argument(
        "--asm_venv", required=True, help="Path to the assembly registry virtual environment."
    )
    args = parser.parse_args()

    try:
        parsed_date = datetime.strptime(args.date, "%m-%d-%Y")
    except ValueError as exc:
        raise ValueError(f"Date '{args.date}' is not in MM-DD-YYYY format") from exc

    date_fmt = datetime.strptime(args.date, "%m-%d-%Y").strftime("%Y-%m-%d")
    gb_registry_flow(
        date=args.date,
        outdir=f"{args.outdir}/{date_fmt}",
        enscode=args.enscode,
        asm_venv=args.asm_venv,
        dry_run=args.dry_run,
    )
