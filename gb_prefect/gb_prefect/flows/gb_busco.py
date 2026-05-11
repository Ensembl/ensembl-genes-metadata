import argparse
import sys
from pathlib import Path
from prefect import flow # type: ignore

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from gb_prefect.tasks.busco import run_nextflow_busco


@flow(name="BUSCO_genome", log_prints=True)
def genome_busco_flow(
    csv_file: str,
    outdir: str,
    dry_run: bool = False,
):
    return run_nextflow_busco(
        csv_file=csv_file,
        outdir=outdir,
        dry_run=dry_run,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true", help="Create the Nextflow command without running it.")
    parser.add_argument("--csv-file", required=True, help="Path to the CSV file.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    args = parser.parse_args()

    genome_busco_flow(csv_file=args.csv_file, 
                      outdir=args.outdir, 
                      dry_run=args.dry_run)
