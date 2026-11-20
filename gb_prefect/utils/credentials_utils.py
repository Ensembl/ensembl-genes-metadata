"""Load the flows' credentials: explicit, or from Prefect Secret blocks.

Flows only need:

    from gb_prefect.utils.credentials_utils import PipelineCredentials, resolve_credentials

    credentials = resolve_credentials(credentials, metadata_secret_block=metadata_secret_block)
    credentials.metadata_params_string   # JSON string, for a pipeline
    credentials.metadata_db_params       # dict, for code connecting to the DB itself
"""

import json
from typing import Any, Optional

from prefect.blocks.system import Secret  # type: ignore

from gb_prefect.models.pipeline_options import PipelineCredentials

__all__ = [
    "DEFAULT_METADATA_SECRET_BLOCK",
    "DEFAULT_SLACK_SECRET_BLOCK",
    "NOT_LOADED_ON_DRY_RUN",
    "PipelineCredentials",
    "resolve_credentials",
]

DEFAULT_METADATA_SECRET_BLOCK = "gb-metadata-db-params"
DEFAULT_SLACK_SECRET_BLOCK = "gb-slack-params"

# Stands in for credentials on a dry run, where nothing is loaded and no params file is written.
NOT_LOADED_ON_DRY_RUN = "<NOT_LOADED_ON_DRY_RUN>"


def _secret_to_json_string(value: Any) -> str:
    """Secret blocks may hold a dict (JSON value) or a plain string; the pipeline wants a JSON string."""
    return value if isinstance(value, str) else json.dumps(value)


def resolve_credentials(
    credentials: Optional[PipelineCredentials] = None,
    *,
    slack_report: bool = False,
    dry_run: bool = False,
    metadata_secret_block: Optional[str] = None,
    slack_secret_block: Optional[str] = None,
) -> PipelineCredentials:
    """Return explicit credentials if given, otherwise load them from Prefect Secret blocks.

    credentials: given explicitly (e.g. when calling a flow by hand); returned as they are.
    slack_report: also load the Slack credentials.
    dry_run: load nothing and return placeholders (no Secret block or server needed).
    metadata_secret_block, slack_secret_block: Secret block names; None means the standard
        ones, DEFAULT_METADATA_SECRET_BLOCK and DEFAULT_SLACK_SECRET_BLOCK.
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
    metadata_block = metadata_secret_block or DEFAULT_METADATA_SECRET_BLOCK
    metadata_params_string = _secret_to_json_string(Secret.load(metadata_block).get())
    slack_params = None
    if slack_report:
        slack_block = slack_secret_block or DEFAULT_SLACK_SECRET_BLOCK
        slack_params = _secret_to_json_string(Secret.load(slack_block).get())
    return PipelineCredentials(metadata_params_string=metadata_params_string, slack_params=slack_params)
