"""`slurm-cli` Prefect worker: submits flow runs as Slurm jobs via `sbatch`.

See gb_prefect/worker/worker.py for the rationale (no Slurm REST API access
on this cluster for this project).
"""

from gb_prefect.worker.worker import SlurmCliWorker

__all__ = ["SlurmCliWorker"]
