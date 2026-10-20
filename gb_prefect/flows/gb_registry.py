import argparse
from datetime import datetime, timedelta
from typing import Optional

from prefect import flow  # type: ignore

from gb_prefect.models.pipeline_options import PipelineCredentials, TaskRunOptions
from gb_prefect.pipelines.assembly_metadata import asm_registry_run
from gb_prefect.tasks.nextflow import run_nextflow_pipeline
from gb_prefect.utils.credentials_utils import (
    DEFAULT_METADATA_SECRET_BLOCK,
    DEFAULT_SLACK_SECRET_BLOCK,
    resolve_credentials,
)
from gb_prefect.utils.enscode_utils import resolve_enscode
from gb_prefect.utils.path_utils import run_folder

# By default, screen NCBI for assemblies released in the last SCREEN_DAYS_BACK days.
SCREEN_DAYS_BACK = 60


@flow(name="gb_registry", log_prints=True)
def gb_registry_flow(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    outdir: str,
    date: Optional[str] = None,
    folder_name: Optional[str] = None,
    enscode: Optional[str] = None,
    dry_run: bool = False,
    credentials: Optional[PipelineCredentials] = None,
    metadata_secret_block: str = DEFAULT_METADATA_SECRET_BLOCK,
):
    """
    Run the assembly registry Nextflow pipeline, screening NCBI for assemblies released
    after date (YYYY-MM-DD; default: SCREEN_DAYS_BACK days before today) and registering the
    ones missing from the metadata DB.
    Output goes to {outdir}/{folder_name}; folder_name defaults to asm_registry_<today>.
    """
    # Work out the date to screen from.
    today = datetime.now()
    if date:
        date = datetime.strptime(date, "%Y-%m-%d").strftime("%Y-%m-%d")
    else:
        date = (today - timedelta(days=SCREEN_DAYS_BACK)).strftime("%Y-%m-%d")
    print(f"Screening NCBI for assemblies released after {date}.")

    # Resolve credentials
    credentials = resolve_credentials(
        credentials, metadata_secret_block, DEFAULT_SLACK_SECRET_BLOCK, slack_report=False, dry_run=dry_run
    )

    # Run the assembly registry pipeline
    run = asm_registry_run(
        outdir=run_folder(outdir, folder_name, f"asm_registry_{today.strftime('%Y-%m-%d')}"),
        enscode=resolve_enscode(enscode, dry_run),
        metadata_params_string=credentials.metadata_params_string,
        date=date,
    )
    return run_nextflow_pipeline(run, TaskRunOptions(dry_run=dry_run))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", required=True, help="Base output directory.")
    parser.add_argument("--enscode", help="Path to the ENSCODE directory; defaults to $ENSCODE.")
    parser.add_argument(
        "--date",
        help="Screen for assemblies released after this date (YYYY-MM-DD). "
        f"Defaults to {SCREEN_DAYS_BACK} days before today.",
    )
    parser.add_argument(
        "--folder-name", help="Output folder name under --outdir. Defaults to asm_registry_<today>."
    )
    parser.add_argument(
        "--metadata-params-string",
        help="JSON string with metadata DB connection parameters. If omitted, credentials are "
        "loaded from the Prefect Secret blocks (see gb_prefect/deployments/create_secrets.py).",
    )
    parser.add_argument(
        "--dry-run", action="store_true", help="Create the Nextflow command without running it."
    )
    args = parser.parse_args()

    gb_registry_flow(
        outdir=args.outdir,
        date=args.date,
        folder_name=args.folder_name,
        enscode=args.enscode,
        dry_run=args.dry_run,
        credentials=(
            PipelineCredentials(metadata_params_string=args.metadata_params_string)
            if args.metadata_params_string
            else None
        ),
    )
