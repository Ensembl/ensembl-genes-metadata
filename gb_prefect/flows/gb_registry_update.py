import argparse
from datetime import datetime
from typing import Optional

from prefect import flow  # type: ignore

from gb_prefect.models.pipeline_options import PipelineCredentials, TaskRunOptions
from gb_prefect.pipelines.assembly_metadata_update import asm_update_run
from gb_prefect.tasks.nextflow import run_nextflow_pipeline
from gb_prefect.utils.credentials_utils import (
    DEFAULT_METADATA_SECRET_BLOCK,
    DEFAULT_SLACK_SECRET_BLOCK,
    resolve_credentials,
)
from gb_prefect.utils.enscode_utils import resolve_enscode
from gb_prefect.utils.path_utils import check_input_file, run_folder


@flow(name="gb_registry_update", log_prints=True)
def gb_registry_update_flow(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    gca_list: str,
    outdir: str,
    slack_report: bool = True,
    folder_name: Optional[str] = None,
    enscode: Optional[str] = None,
    dry_run: bool = False,
    credentials: Optional[PipelineCredentials] = None,
    metadata_secret_block: str = DEFAULT_METADATA_SECRET_BLOCK,
    slack_secret_block: str = DEFAULT_SLACK_SECRET_BLOCK,
):
    """Run the assembly metadata update Nextflow pipeline for the GCA list (the whole file, as
    given). See gb_registry_update_by_date_flow for screening by release date instead.

    gca_list must be an absolute path. Output goes to {outdir}/{folder_name}; folder_name
    defaults to asm_update_<today>.
    """
    # Resolve credentials
    credentials = resolve_credentials(
        credentials, metadata_secret_block, slack_secret_block, slack_report, dry_run=dry_run
    )
    
    # Run the assembly metadata update pipeline for the specified GCA accessions.
    run = asm_update_run(
        outdir=run_folder(outdir, folder_name, f"asm_update_{datetime.now().strftime('%Y-%m-%d')}"),
        enscode=resolve_enscode(enscode, dry_run),
        metadata_params_string=credentials.metadata_params_string,
        gca_list=check_input_file(gca_list, dry_run),
        slack_report=slack_report,
        slack_params=credentials.slack_params,
    )
    return run_nextflow_pipeline(run, TaskRunOptions(dry_run=dry_run))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--gca-list", required=True, help="Absolute path to the file listing the GCA accessions to update."
    )
    parser.add_argument("--outdir", required=True, help="Base output directory.")
    parser.add_argument("--enscode", help="Path to the ENSCODE directory; defaults to $ENSCODE.")
    parser.add_argument(
        "--folder-name", help="Output folder name under --outdir. Defaults to asm_update_<today>."
    )
    parser.add_argument(
        "--slack-report", action="store_true", help="Send a Slack report after the Nextflow run."
    )
    parser.add_argument(
        "--slack-params",
        help="JSON string with Slack bot connection parameters, used with --slack-report.",
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

    gb_registry_update_flow(
        gca_list=args.gca_list,
        outdir=args.outdir,
        slack_report=args.slack_report,
        folder_name=args.folder_name,
        enscode=args.enscode,
        dry_run=args.dry_run,
        credentials=(
            PipelineCredentials(
                metadata_params_string=args.metadata_params_string,
                slack_params=args.slack_params,
            )
            if args.metadata_params_string
            else None
        ),
    )
