"""Facts about the current flow run that stay the same when it is retried.

Defaults derived from "now" (a screen date, an output folder) change if a retry starts on a
later day, so a retry would look at different inputs, or write to a new folder where Nextflow
has nothing to -resume. These come from the flow run itself instead.
"""

from datetime import datetime

from prefect.runtime import flow_run  # type: ignore


def run_start() -> datetime:
    """The flow run's expected start time, in local time (naive). Prefect sets it once, when
    the run is first scheduled, so it doesn't change on a retry. Outside a flow run: now."""
    return datetime.fromtimestamp(flow_run.scheduled_start_time.timestamp())


def run_name_suffix() -> str:
    """"_<flow run name>" (e.g. "_sigma-normae"), unique per run and kept on a retry; "" outside
    a flow run."""
    name = flow_run.name
    return f"_{name}" if name else ""


def is_retry() -> bool:
    """True if this flow run has run before (a retry, automatic or from the UI)."""
    return (flow_run.run_count or 0) > 1
