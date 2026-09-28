"""Minimal flow to validate the `slurm-cli` worker end to end.

Deploy this to a real `slurm-cli` work pool and trigger one run to confirm:
the worker actually submits a Slurm job via `sbatch`, the job's environment
has `gb_prefect`/`prefect` importable (via the deployment's `setup_commands`),
and the flow run reports its own completion back to Prefect normally.

Not a production flow -- delete or move once the worker is proven.
"""

import os
import socket

from prefect import flow, task  # type: ignore


@task
def report_environment():
    return {
        "hostname": socket.gethostname(),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "slurm_job_partition": os.environ.get("SLURM_JOB_PARTITION"),
        "cwd": os.getcwd(),
    }


@flow(name="slurm_cli_smoke_test", log_prints=True)
def slurm_cli_smoke_test_flow():
    """Prints the environment it actually ran in, to confirm it executed inside a real Slurm job."""
    info = report_environment()
    print(f"Running on host: {info['hostname']}")
    print(f"SLURM_JOB_ID: {info['slurm_job_id']}")
    print(f"SLURM_JOB_PARTITION: {info['slurm_job_partition']}")
    print(f"Working directory: {info['cwd']}")

    if info["slurm_job_id"] is None:
        raise RuntimeError(
            "SLURM_JOB_ID is not set -- this flow run did not execute inside a Slurm job. "
            "Check the pool's type is 'slurm-cli' and the worker submitting it is running "
            "on a host with sbatch access."
        )

    return info


if __name__ == "__main__":
    slurm_cli_smoke_test_flow()
