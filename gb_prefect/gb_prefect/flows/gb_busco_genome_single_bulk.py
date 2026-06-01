import argparse
from datetime import datetime
import pathlib
import sys
from prefect import flow # type: ignore
from typing import Optional

PACKAGE_ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from gb_prefect.tasks.busco import run_nextflow_busco
from gb_prefect.utils.stats_split_csv import split_csv

@flow(name="BUSCO_genome_single_bulk", log_prints=True)
def genome_busco_master_flow(
    csv_file: str,
    outdir: str,
    asm_venv: str,
    enscode: Optional[str] = None,
    create_artifact: bool = True,
    dry_run: bool = False):

    nxf_outdir = str(pathlib.Path(outdir) / f"busco_genome_{datetime.now().strftime('%Y%m%d_%H%M%S')}")

    csv_files = split_csv(csv_file, outdir)
    futures = [
        run_nextflow_busco.submit(
            csv_file=f,
            outdir=nxf_outdir,
            asm_venv=asm_venv,
            enscode=enscode,
            create_artifact=create_artifact,
            dry_run=dry_run,
        )
        for f in csv_files
    ]
    return [f.result() for f in futures]

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv-file", required=True, help="Path to the CSV file.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--asm-venv", required=True, help="Path to the ASM virtual environment.")
    parser.add_argument("--enscode", help="ENSCODE value.")
    parser.add_argument("--create-artifact", action="store_true", help="Create artifact after running the pipeline.")
    parser.add_argument("--dry-run", action="store_true", help="Create the Nextflow command without running it.")
    args = parser.parse_args()

    genome_busco_master_flow(
        csv_file=args.csv_file,
        outdir=args.outdir,
        asm_venv=args.asm_venv,
        enscode=args.enscode,
        create_artifact=args.create_artifact,
        dry_run=args.dry_run,
    )
