"""Shared parameter bundles for the assembly registry/update Prefect tasks and flows."""

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class PipelineCredentials:
    """Credential JSON strings passed through to a Nextflow pipeline run.

    Attributes:
        metadata_params_string: JSON string with metadata database connection parameters.
        slack_params: JSON string with Slack bot connection parameters. Required only when
            Slack reporting is enabled.
    """

    metadata_params_string: str
    slack_params: Optional[str] = None


@dataclass(frozen=True)
class TaskRunOptions:
    """Execution-control flags for a Prefect task, separate from its pipeline inputs.

    Attributes:
        dry_run: If True, build the sbatch script without submitting it.
        create_artifact: If True, create a Prefect artifact summarizing the run.
    """

    dry_run: bool = False
    create_artifact: bool = True
