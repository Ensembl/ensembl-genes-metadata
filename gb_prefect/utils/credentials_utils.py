import json
from typing import Any, Dict, Optional

from gb_prefect.models.pipeline_options import PipelineCredentials

DEFAULT_METADATA_SECRET_BLOCK = "gb-metadata-db-params"
DEFAULT_SLACK_SECRET_BLOCK = "gb-slack-params"

REDACTED = "<REDACTED>"


def _secret_to_json_string(value: Any) -> str:
    """Secret blocks may hold a dict (JSON value) or a plain string; the pipeline wants a JSON string."""
    return value if isinstance(value, str) else json.dumps(value)


def resolve_credentials(
    credentials: Optional[PipelineCredentials],
    metadata_secret_block: str,
    slack_secret_block: str,
    slack_report: bool,
) -> PipelineCredentials:
    """Return explicit credentials if given, otherwise load them from Prefect Secret blocks.

    Explicit credentials are used by the standalone CLI; deployments leave them unset so the
    values live in Secret blocks instead of in the deployment's (plain-text) parameters.
    """
    if credentials is not None:
        return credentials

    from prefect.blocks.system import Secret  # type: ignore  # pylint: disable=import-outside-toplevel

    metadata_params_string = _secret_to_json_string(Secret.load(metadata_secret_block).get())
    slack_params = None
    if slack_report:
        slack_params = _secret_to_json_string(Secret.load(slack_secret_block).get())
    return PipelineCredentials(metadata_params_string=metadata_params_string, slack_params=slack_params)


def load_metadata_db_params(
    credentials: Optional[PipelineCredentials],
    metadata_secret_block: str = DEFAULT_METADATA_SECRET_BLOCK,
) -> Dict[str, Any]:
    """Return the assembly metadata DB connection params as a dict, from credentials if given,
    otherwise from the metadata_secret_block Prefect Secret block."""
    credentials = resolve_credentials(
        credentials, metadata_secret_block, DEFAULT_SLACK_SECRET_BLOCK, slack_report=False
    )
    return json.loads(credentials.metadata_params_string)


def redact_credentials(text: str, credentials: PipelineCredentials) -> str:
    """Replace credential JSON strings in text, so it can be logged or published as an artifact."""
    for secret in (credentials.metadata_params_string, credentials.slack_params):
        if secret:
            text = text.replace(secret, REDACTED)
    return text
