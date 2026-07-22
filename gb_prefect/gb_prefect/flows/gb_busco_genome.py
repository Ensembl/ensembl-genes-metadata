import argparse
import sys
from pathlib import Path
from prefect import flow # type: ignore
from typing import Optional

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from gb_prefect.tasks.busco import run_nextflow_busco


@flow(name="BUSCO_genome", log_prints=True)
def genome_busco_flow(
    csv_file: str,
    outdir: str,
    asm_venv: str,
    enscode: Optional[str] = None,
    dry_run: bool = False,
    create_artifact: bool = True,   
):
    return run_nextflow_busco(
        csv_file=csv_file,
        outdir=outdir,
        asm_venv=asm_venv,
        enscode=enscode,
        dry_run=dry_run,
        create_artifact=create_artifact,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv-file", required=True, help="Path to the CSV file.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--asm-venv", required=True, help="Path to the asm_metadata virtual environment.")
    parser.add_argument("--enscode", help="ENSCODE value.")
    parser.add_argument("--dry-run", action="store_true", help="Create the Nextflow command without running it.")
    parser.add_argument("--create-artifact", action="store_true", help="Create artifact after running the pipeline.")
    args = parser.parse_args()

    genome_busco_flow(csv_file=args.csv_file, 
                      outdir=args.outdir, 
                      asm_venv=args.asm_venv,
                      enscode=args.enscode,
                      dry_run=args.dry_run)
