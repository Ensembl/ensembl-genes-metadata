import argparse
from datetime import datetime
from typing import Optional

from prefect import flow  # type: ignore
from gb_prefect.models.pipeline_options import PipelineCredentials, TaskRunOptions
from gb_prefect.utils.credentials_utils import (
    DEFAULT_METADATA_SECRET_BLOCK,
    DEFAULT_SLACK_SECRET_BLOCK,
    resolve_credentials,
)
from gb_prefect.tasks.registry_update import update_assemblies


@flow(name="gb_registry_update", log_prints=True)
def gb_registry_update_flow(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    gca_list: str,
    outdir: str,
    asm_venv: str,
    slack_report: bool = True,
    date: Optional[str] = None,
    enscode: Optional[str] = None,
    dry_run: bool = False,
    credentials: Optional[PipelineCredentials] = None,
    metadata_secret_block: str = DEFAULT_METADATA_SECRET_BLOCK,
    slack_secret_block: str = DEFAULT_SLACK_SECRET_BLOCK,
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
        credentials=resolve_credentials(
            credentials, metadata_secret_block, slack_secret_block, slack_report
        ),
        slack_report=slack_report,
        date=date,
        enscode=enscode,
        run_options=TaskRunOptions(dry_run=dry_run),
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
        required=False,
        help="JSON string with metadata database connection parameters. If omitted, credentials are "
        "loaded from the Prefect Secret blocks (see gb_prefect/deployments/create_secrets.py).",
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
        outdir=args.outdir,
        gca_list=args.gca_list,
        enscode=args.enscode,
        asm_venv=args.asm_venv,
        credentials=(
            PipelineCredentials(
                metadata_params_string=args.metadata_params_string,
                slack_params=args.slack_params,
            )
            if args.metadata_params_string
            else None
        ),
        slack_report=args.slack_report,
        dry_run=args.dry_run,
    )
