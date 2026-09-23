"""Logging filter that redacts Prefect's own secrets from worker logs.

This protects credentials that flow through `SlurmCliWorkerConfiguration.env`
into the generated sbatch script and into the worker's own log output --
namely Prefect's API auth (PREFECT_API_KEY, PREFECT_API_AUTH_STRING), which
every flow run needs to reach the Prefect API.

It does NOT redact database/Slack credentials that individual tasks embed in
their own Nextflow commands (see gb_prefect/tasks/registry_update.py) -- that
is a separate, task-level concern this worker package does not touch.
"""

import logging
import os
from typing import Any, Iterable, Mapping

from prefect.logging.filters import redact_substr

from gb_prefect.worker.config import LOG_MASK_PATTERNS


def mask_sensitive_data(
    value: Any, env: Mapping[str, Any], patterns: Iterable[str] = LOG_MASK_PATTERNS
) -> Any:
    secrets = (env.get(pattern) or pattern for pattern in patterns)

    for secret in sorted(filter(None, secrets), key=len, reverse=True):
        value = redact_substr(value, secret)
    return value


class RedactingFilter(logging.Filter):
    def filter(self, record: logging.LogRecord) -> bool:
        record.msg = mask_sensitive_data(record.msg, env=os.environ)
        if record.args:
            record.args = tuple(
                mask_sensitive_data(arg, os.environ) for arg in record.args
            )

        return True
