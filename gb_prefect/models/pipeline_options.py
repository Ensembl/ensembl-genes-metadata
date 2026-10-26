"""Shared parameter bundles for the Prefect tasks and flows."""

import json
from dataclasses import dataclass
from typing import Any, Dict, Optional


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

    @property
    def metadata_db_params(self) -> Dict[str, Any]:
        """The metadata DB connection parameters as a dict, for code connecting to the DB itself
        (e.g. with pymysql); pipelines get metadata_params_string as it is."""
        return json.loads(self.metadata_params_string)


@dataclass(frozen=True)
class TaskRunOptions:
    """Execution-control flags for a Prefect task.

    Attributes:
        dry_run: If True, build the sbatch script without submitting it.
        create_artifact: If True, create a Prefect artifact summarizing the run.
    """

    dry_run: bool = False
    create_artifact: bool = True
