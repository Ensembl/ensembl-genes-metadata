import argparse
from datetime import datetime
from typing import Optional

from prefect import flow  # type: ignore
from gb_prefect.tasks.registry_update import update_assemblies


@flow(name="gb_registry_update", log_prints=True)
def gb_registry_update_flow(
    gca_list: str,
    outdir: str,
    asm_venv: str,
    metadata_params_string: str,
    slack_report: bool = True,
    slack_params: Optional[str] = None,
    date: Optional[str] = None,
    enscode: Optional[str] = None,
    dry_run: bool = False,
):
    """Run the assembly metadata update Nextflow pipeline for the given GCA list."""
    if not date:
        date = datetime.now().strftime("%Y-%m-%d")
    else:
        date = datetime.strptime(date, "%Y-%m-%d").strftime("%Y-%m-%d")

    return update_assemblies(
        gca_list=gca_list,
        outdir=f"{outdir}/asm_update_{date}",
        asm_venv=asm_venv,
        metadata_params_string=metadata_params_string,
        slack_report=slack_report,
        slack_params=slack_params,
        date=date,
        enscode=enscode,
        dry_run=dry_run,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--gca-list", required=True, help="Path to the file containing the list of GCA accessions to update."
    )
    parser.add_argument("--outdir", required=True, help="Base output directory.")
    parser.add_argument(
        "--asm_venv", required=True, help="Path to the assembly registry virtual environment."
    )
    parser.add_argument(
        "--slack-report", action="store_true", help="Whether to send a Slack report after the Nextflow run."
    )
    parser.add_argument(
        "--slack-params",
        required=False,
        help="JSON string with Slack bot connection parameters. Required if --slack-report is set.",
    )
    parser.add_argument(
        "--metadata-params-string",
        required=True,
        help="JSON string with metadata database connection parameters.",
    )
    parser.add_argument("--enscode", required=True, help="Path to ENSCODE directory.")
    parser.add_argument("--date", required=False, help="Date for the registry run (e.g., YYYY-MM-DD).")
    parser.add_argument(
        "--dry-run", action="store_true", help="Create the Nextflow command without running it."
    )
    args = parser.parse_args()

    if not args.date:
        run_date = datetime.now().strftime("%Y-%m-%d")
    else:
        run_date = datetime.strptime(args.date, "%Y-%m-%d").strftime("%Y-%m-%d")

    gb_registry_update_flow(
        date=run_date,
        outdir=f"{args.outdir}/asm_update_{run_date}",
        gca_list=args.gca_list,
        enscode=args.enscode,
        asm_venv=args.asm_venv,
        metadata_params_string=args.metadata_params_string,
        slack_report=args.slack_report,
        slack_params=args.slack_params,
        dry_run=args.dry_run,
    )
