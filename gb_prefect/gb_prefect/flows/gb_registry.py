import argparse
import sys
from pathlib import Path
from prefect import flow # type: ignore
from datetime import datetime


PACKAGE_ROOT = Path(__file__).resolve().parents[2]
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from gb_prefect.tasks.registry import register_assemblies

@flow(name="gb_registry", log_prints=True)
def gb_registry_flow(
    date: str,
    outdir: str,
    enscode: str,
    dry_run: bool = False,
):

    return register_assemblies(
        date=date,
        outdir=f"{outdir}/{date}",
        enscode=enscode,
        dry_run=dry_run,
    )

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true", help="Create the Nextflow command without running it.")
    parser.add_argument("--date", required=True, help="Date for the registry run (e.g., MM-DD-YYYY).")
    parser.add_argument("--outdir", required=True, help="Base output directory.")
    parser.add_argument("--enscode", required=True, help="Path to ENSCODE directory.")
    args = parser.parse_args()

    try:
        parsed_date = datetime.strptime(args.date, "%m-%d-%Y")
    except ValueError:
        raise ValueError(f"Date '{args.date}' is not in MM-DD-YYYY format")

    date_fmt = datetime.strptime(args.date, "%m-%d-%Y").strftime("%Y-%m-%d")
    gb_registry_flow(date=args.date,
                     outdir=f"{args.outdir}/{date_fmt}",
                     enscode=args.enscode,
                     dry_run=args.dry_run)