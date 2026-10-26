import json
from typing import Any, Optional
from prefect.blocks.system import Secret  # type: ignore  # pylint: disable=import-outside-toplevel
from gb_prefect.models.pipeline_options import PipelineCredentials

DEFAULT_METADATA_SECRET_BLOCK = "gb-metadata-db-params"
DEFAULT_SLACK_SECRET_BLOCK = "gb-slack-params"

# Stands in for credentials on a dry run, where nothing is loaded and no params file is written.
NOT_LOADED_ON_DRY_RUN = "<NOT_LOADED_ON_DRY_RUN>"


def _secret_to_json_string(value: Any) -> str:
    """Secret blocks may hold a dict (JSON value) or a plain string; the pipeline wants a JSON string."""
    return value if isinstance(value, str) else json.dumps(value)


def resolve_credentials(
    credentials: Optional[PipelineCredentials],
    metadata_secret_block: str,
    slack_secret_block: str,
    slack_report: bool,
    dry_run: bool = False,
) -> PipelineCredentials:
    """Return explicit credentials if given, otherwise load them from Prefect Secret blocks.
    dry_run=True loads nothing and returns placeholders.
    """
    # Return explicit credentials if given
    if credentials is not None:
        return credentials

    # Return placeholders if dry_run
    if dry_run:
        return PipelineCredentials(
            metadata_params_string=NOT_LOADED_ON_DRY_RUN,
            slack_params=NOT_LOADED_ON_DRY_RUN if slack_report else None,
        )

    # Return credentials from secret blocks
    metadata_params_string = _secret_to_json_string(Secret.load(metadata_secret_block).get())
    slack_params = None
    if slack_report:
        slack_params = _secret_to_json_string(Secret.load(slack_secret_block).get())
    return PipelineCredentials(metadata_params_string=metadata_params_string, slack_params=slack_params)
