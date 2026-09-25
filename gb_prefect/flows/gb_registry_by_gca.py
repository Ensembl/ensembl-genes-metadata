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
from gb_prefect.tasks.registry import register_assemblies


@flow(name="gb_registry_by_gca", log_prints=True)
def gb_registry_by_gca_flow(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    gca_list: str,
    outdir: str,
    asm_venv: str,
    enscode: Optional[str] = None,
    dry_run: bool = False,
    credentials: Optional[PipelineCredentials] = None,
    metadata_secret_block: str = DEFAULT_METADATA_SECRET_BLOCK,
):
    """Run the assembly registry Nextflow pipeline for the GCA accessions listed in gca_list
    (--add_gca mode), instead of screening NCBI (see gb_registry_flow for that mode)."""
    return register_assemblies(
        outdir=f"{outdir}/asm_registry_gca_{datetime.now().strftime('%Y-%m-%d')}",
        asm_venv=asm_venv,
        credentials=resolve_credentials(
            credentials, metadata_secret_block, DEFAULT_SLACK_SECRET_BLOCK, slack_report=False
        ),
        gca_list=gca_list,
        enscode=enscode,
        run_options=TaskRunOptions(dry_run=dry_run),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--gca-list",
        required=True,
        help="Path to a file with the GCA accessions to register (lines starting with GCA_).",
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
    parser.add_argument(
        "--dry-run", action="store_true", help="Create the Nextflow command without running it."
    )
    args = parser.parse_args()

    gb_registry_by_gca_flow(
        gca_list=args.gca_list,
        outdir=args.outdir,
        asm_venv=args.asm_venv,
        enscode=args.enscode,
        dry_run=args.dry_run,
        credentials=(
            PipelineCredentials(metadata_params_string=args.metadata_params_string)
            if args.metadata_params_string
            else None
        ),
    )
