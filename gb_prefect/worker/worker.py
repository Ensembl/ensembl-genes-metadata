"""Prefect worker that submits flow runs as Slurm jobs via the `sbatch` CLI.

This cluster has `slurmrestd`/JWT token issuance unavailable for this project
(the Codon Submitter VM was specifically provisioned for direct Slurm CLI
access instead), so this does not use Slurm's REST API the way prefect-slurm's
SlurmWorker does. Job submission goes through `sbatch` directly, run as an
async subprocess so the worker's polling loops are never blocked waiting on
a Slurm job to finish.

As with any Prefect worker, the whole flow run -- including any task that
shells out to Nextflow -- executes inside a single Slurm job; tasks are not
separate Slurm jobs.
"""

import logging
import os
import re
from typing import Optional

import anyio
from anyio.abc import TaskStatus
from prefect.client.schemas import FlowRun
from prefect.exceptions import InfrastructureError
from prefect.logging.loggers import PrefectLogAdapter
from prefect.workers.base import BaseWorker, BaseWorkerResult
from tenacity import retry, stop_after_attempt, wait_fixed, wait_random

from gb_prefect.worker.config import (
    SlurmCliWorkerConfiguration,
    SlurmCliWorkerTemplateVariables,
)
from gb_prefect.worker.log_filters import RedactingFilter

MAX_ATTEMPTS = int(os.getenv("GB_PREFECT_WORKER_MAX_ATTEMPTS", 3))
RETRY_MIN_DELAY_SECONDS = int(os.getenv("GB_PREFECT_WORKER_RETRY_MIN_DELAY_SECONDS", 10))
RETRY_MIN_DELAY_JITTER_SECONDS = int(
    os.getenv("GB_PREFECT_WORKER_RETRY_MIN_DELAY_JITTER_SECONDS", 0)
)
RETRY_MAX_DELAY_JITTER_SECONDS = int(
    os.getenv("GB_PREFECT_WORKER_RETRY_MAX_DELAY_JITTER_SECONDS", 20)
)

SLURM_JOB_ID_PATTERN = re.compile(r"Submitted batch job (\d+)")


class SlurmCliWorker(
    BaseWorker[SlurmCliWorkerConfiguration, SlurmCliWorkerTemplateVariables, BaseWorkerResult]
):
    """A Prefect worker that submits flow runs as Slurm jobs via `sbatch`."""

    type: str = "slurm-cli"
    job_configuration = SlurmCliWorkerConfiguration
    job_configuration_variables = SlurmCliWorkerTemplateVariables
    _documentation_url = (
        "https://github.com/Ensembl/ensembl-genes-metadata/tree/main/gb_prefect"
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Attach filter to logger and guard against attaching same filter twice to reusable loggers
        if isinstance(self._logger, logging.Logger):
            actual_logger = self._logger
        else:
            actual_logger = self._logger.logger

        if not any(isinstance(f, RedactingFilter) for f in actual_logger.filters):
            actual_logger.addFilter(RedactingFilter())

    def get_flow_run_logger(self, flow_run: FlowRun) -> PrefectLogAdapter:
        adapter = super().get_flow_run_logger(flow_run)

        if not any(isinstance(f, RedactingFilter) for f in adapter.logger.filters):
            adapter.logger.addFilter(RedactingFilter())

        return adapter

    async def run(
        self,
        flow_run: FlowRun,
        configuration: SlurmCliWorkerConfiguration,
        task_status: Optional[TaskStatus[str]] = None,
    ) -> BaseWorkerResult:
        logger = self.get_flow_run_logger(flow_run)

        configuration.working_dir.mkdir(parents=True, exist_ok=True)
        script_path = configuration.working_dir / f"slurm_cli_{flow_run.id}.sh"
        # Secrets (e.g. PREFECT_API_KEY) are exported in this script, so keep it
        # readable only by the submitting user.
        script_path.write_text(configuration.script)
        script_path.chmod(0o700)

        logger.info(f"Submitting flow run {flow_run.id} to Slurm via sbatch")
        job_id = await self._submit_sbatch(script_path, configuration.working_dir)
        logger.info(f"Submitted flow run {flow_run.id} to Slurm job {job_id}")

        if task_status:
            task_status.started(job_id)

        return BaseWorkerResult(status_code=0, identifier=job_id)

    @retry(
        stop=stop_after_attempt(MAX_ATTEMPTS),
        wait=wait_fixed(RETRY_MIN_DELAY_SECONDS)
        + wait_random(
            RETRY_MIN_DELAY_JITTER_SECONDS,
            RETRY_MAX_DELAY_JITTER_SECONDS,
        ),
        reraise=True,
    )
    async def _submit_sbatch(self, script_path, working_dir) -> str:
        result = await anyio.run_process(
            ["sbatch", str(script_path)], cwd=str(working_dir), check=False
        )

        if result.returncode != 0:
            raise InfrastructureError(
                f"sbatch submission failed (exit {result.returncode}): "
                f"{result.stderr.decode(errors='replace')}"
            )

        stdout = result.stdout.decode(errors="replace")
        match = SLURM_JOB_ID_PATTERN.search(stdout)
        if not match:
            raise InfrastructureError(
                f"Could not parse Slurm job ID from sbatch output: {stdout!r}"
            )

        return match.group(1)
