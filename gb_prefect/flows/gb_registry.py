import argparse
from datetime import datetime, timedelta
from typing import Optional

from prefect import flow  # type: ignore
from gb_prefect.models.pipeline_options import PipelineCredentials, TaskRunOptions
from gb_prefect.utils.credentials_utils import (
    DEFAULT_METADATA_SECRET_BLOCK,
    DEFAULT_SLACK_SECRET_BLOCK,
    resolve_credentials,
)
from gb_prefect.tasks.registry import register_assemblies

# By default, screen NCBI for assemblies released in the last SCREEN_DAYS_BACK days.
SCREEN_DAYS_BACK = 60


@flow(name="gb_registry", log_prints=True)
def gb_registry_flow(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    outdir: str,
    asm_venv: str,
    date: Optional[str] = None,
    enscode: Optional[str] = None,
    dry_run: bool = False,
    credentials: Optional[PipelineCredentials] = None,
    metadata_secret_block: str = DEFAULT_METADATA_SECRET_BLOCK,
):
    """Run the assembly registry Nextflow pipeline, screening NCBI for assemblies released
    after date (MM-DD-YYYY; default: SCREEN_DAYS_BACK days before today) and registering the
    ones missing from the metadata DB. Output goes to a folder named after today's date.
    """
    today = datetime.now()
    if date:
        date = datetime.strptime(date, "%m-%d-%Y").strftime("%m-%d-%Y")
    else:
        date = (today - timedelta(days=SCREEN_DAYS_BACK)).strftime("%m-%d-%Y")
    print(f"Screening NCBI for assemblies released after {date}.")

    return register_assemblies(
        outdir=f"{outdir}/asm_registry_{today.strftime('%Y-%m-%d')}",
        asm_venv=asm_venv,
        credentials=resolve_credentials(
            credentials, metadata_secret_block, DEFAULT_SLACK_SECRET_BLOCK, slack_report=False
        ),
        date=date,
        enscode=enscode,
        run_options=TaskRunOptions(dry_run=dry_run),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dry-run", action="store_true", help="Create the Nextflow command without running it."
    )
    parser.add_argument(
        "--date",
        required=False,
        help="Screen for assemblies released after this date (MM-DD-YYYY). "
        f"Defaults to {SCREEN_DAYS_BACK} days before today.",
    )
    parser.add_argument("--outdir", required=True, help="Base output directory.")
    parser.add_argument("--enscode", required=True, help="Path to ENSCODE directory.")
    parser.add_argument(
        "--asm_venv", required=True, help="Path to the assembly registry virtual environment."
    )
    parser.add_argument(
        "--metadata-params-string",
        required=False,
        help="JSON string with metadata database connection parameters. If omitted, credentials are "
        "loaded from the Prefect Secret blocks (see gb_prefect/deployments/create_secrets.py).",
    )
    args = parser.parse_args()

    gb_registry_flow(
        outdir=args.outdir,
        asm_venv=args.asm_venv,
        date=args.date,
        enscode=args.enscode,
        dry_run=args.dry_run,
        credentials=(
            PipelineCredentials(metadata_params_string=args.metadata_params_string)
            if args.metadata_params_string
            else None
        ),
    )
