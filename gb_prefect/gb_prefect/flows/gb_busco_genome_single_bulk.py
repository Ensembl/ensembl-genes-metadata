import argparse
import pathlib
import sys  
from prefect import flow # type: ignore

PACKAGE_ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from gb_prefect.flows.gb_busco_genome import genome_busco_flow
from gb_prefect.utils.stats_split_csv import split_csv

@flow(name="BUSCO_genome_single_bulk", log_prints=True)
def genome_busco_master_flow(csv_file: str, outdir: str, dry_run: bool = False):
    csv_files = split_csv(csv_file, outdir)
    futures = [
        genome_busco_flow(csv_file=f, outdir=outdir, dry_run=dry_run, return_state=False)
        for f in csv_files       
    ]
    return futures

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true", help="Create the Nextflow command without running it.")
    parser.add_argument("--csv-file", required=True, help="Path to the CSV file.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    args = parser.parse_args()

    genome_busco_master_flow(csv_file=args.csv_file, 
                             outdir=args.outdir, 
                             dry_run=args.dry_run)