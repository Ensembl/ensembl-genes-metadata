import argparse
import sys
from pathlib import Path
from prefect import flow # type: ignore
from datetime import datetime


PACKAGE_ROOT = Path(__file__).resolve().parents[2]
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from gb_prefect.tasks.registry_update import update_assemblies

@flow(name="gb_registry_update", log_prints=True)
def gb_registry_update_flow(
    gca_list: str,
    outdir: str,
    asm_venv: str,
    slack_report: bool = True,
    date: str = None,
    enscode: str = None,
    dry_run: bool = False,
):
    
    return update_assemblies(
        gca_list=gca_list,
        outdir=f"{outdir}/{date}",
        asm_venv=asm_venv,
        slack_report=slack_report,
        date=date,
        enscode=enscode,
        dry_run=dry_run,
    )

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--gca-list", required=True, help="Path to the file containing the list of GCA accessions to update.")
    parser.add_argument("--outdir", required=True, help="Base output directory.")
    parser.add_argument("--asm_venv", required=True, help="Path to the assembly registry virtual environment.")
    parser.add_argument("--slack-report", action="store_true", help="Whether to send a Slack report after the Nextflow run.")
    parser.add_argument("--enscode", required=True, help="Path to ENSCODE directory.")
    parser.add_argument("--date", required=True, help="Date for the registry run (e.g., MM-DD-YYYY).")
    parser.add_argument("--dry-run", action="store_true", help="Create the Nextflow command without running it.")
    args = parser.parse_args()
    
    if not args.date:
        date = datetime.now().strftime("%Y-%m-%d")

    gb_registry_update_flow(date=args.date,
                     outdir=f"{args.outdir}/asm_update_{date}",
                     gca_list=args.gca_list,
                     enscode=args.enscode,
                     asm_venv=args.asm_venv,
                     slack_report=args.slack_report,
                     dry_run=args.dry_run)