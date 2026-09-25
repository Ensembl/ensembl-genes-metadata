"""Job configuration for the `slurm-cli` Prefect worker.

Builds one sbatch script per flow run and hands it to the worker to submit via
the `sbatch` CLI. Unlike prefect-slurm's SlurmWorker (which submits through
Slurm's REST API), there is no REST payload here -- resource requests are
rendered directly as `#SBATCH` directives in the generated script, since that
is what `sbatch <script>` actually consumes.
"""

import re
import shlex
from pathlib import Path
from typing import TYPE_CHECKING, List, Optional
from uuid import UUID

from prefect.client.schemas import FlowRun
from prefect.workers.base import BaseJobConfiguration, BaseVariables
from pydantic import Field

if TYPE_CHECKING:
    from prefect.client.schemas.objects import WorkPool
    from prefect.client.schemas.responses import DeploymentResponse
    from prefect.flows import Flow as APIFlow


LOG_MASK_PATTERNS = (
    "PREFECT_API_AUTH_STRING",
    "PREFECT_API_KEY",
)

_JOB_NAME_SANITIZE_PATTERN = re.compile(r"[^A-Za-z0-9_.-]")


class SlurmCliWorkerConfiguration(BaseJobConfiguration):
    """Configuration for the `slurm-cli` worker.

    The whole flow run (including any task that shells out to Nextflow) is
    submitted as a single Slurm job -- tasks are not separate Slurm jobs.
    """

    cpu: int = Field(default=1, description="CPU count required for the flow")
    memory: int = Field(default=4, description="Memory in GB required for the flow")
    partition: Optional[str] = Field(default=None, description="Slurm partition to use")
    account: Optional[str] = Field(default=None, description="Slurm account to charge the job to")
    qos: Optional[str] = Field(default=None, description="Slurm QOS to submit the job under")
    time_limit: int = Field(
        default=1,
        title="Time Limit",
        description="Max wall time in hours for the flow",
    )
    working_dir: Path = Field(
        title="Working Directory",
        description="Directory to run the flow from, and where the sbatch script and Slurm "
        "output/error files are written",
    )
    setup_commands: List[str] = Field(
        default_factory=list,
        title="Setup Commands",
        description=(
            "Shell commands run before the flow-run command, e.g. "
            "'module load nextflow/24.10.3', 'source /path/to/asm_venv/bin/activate'. "
            "The resulting environment must have this package and prefect installed, "
            "since the flow run itself executes inside this Slurm job."
        ),
        examples=[["module load nextflow/24.10.3", "source ~/asm_venv/bin/activate"]],
    )
    shebang: str = Field(
        default="#!/bin/bash",
        pattern=r"^#!/.+$",
        description="Indicates which shell to use when running the slurm job",
    )
    script: str = Field(
        default="",
        description="Generated sbatch script for the flow run",
        json_schema_extra=dict(template="{{command}}"),
    )

    def prepare_for_flow_run(
        self,
        flow_run: FlowRun,
        deployment: Optional["DeploymentResponse"] = None,
        flow: Optional["APIFlow"] = None,
        work_pool: Optional["WorkPool"] = None,
        worker_name: Optional[str] = None,
        worker_id: Optional[UUID] = None,
    ):
        """Builds the sbatch script (self.script) to submit for this flow run."""
        super().prepare_for_flow_run(
            flow_run=flow_run,
            deployment=deployment,
            flow=flow,
            work_pool=work_pool,
            worker_name=worker_name,
            worker_id=worker_id,
        )

        script_segments = [
            self.shebang.strip(),
            self._sbatch_directives_segment(flow_run),
            self._env_export_segment(),
            self._setup_commands_segment(),
            self.command,
        ]
        self.script = "\n".join(segment for segment in script_segments if segment)

    def _sbatch_directives_segment(self, flow_run: FlowRun) -> str:
        lines = [
            f"#SBATCH --job-name={self._sanitized_job_name(flow_run)}",
            f"#SBATCH --chdir={self.working_dir}",
            f"#SBATCH --output={self.working_dir}/slurm_%j.out",
            f"#SBATCH --error={self.working_dir}/slurm_%j.err",
            f"#SBATCH --cpus-per-task={self.cpu}",
            f"#SBATCH --mem={self.memory}G",
            f"#SBATCH --time={self.time_limit:02d}:00:00",
        ]
        if self.partition:
            lines.append(f"#SBATCH --partition={self.partition}")
        if self.account:
            lines.append(f"#SBATCH --account={self.account}")
        if self.qos:
            lines.append(f"#SBATCH --qos={self.qos}")

        return "\n".join(lines)

    def _env_export_segment(self) -> Optional[str]:
        if not self.env:
            return None

        return "\n".join(
            f"export {key}={shlex.quote(str(value))}"
            for key, value in self.env.items()
            if value is not None
        )

    def _setup_commands_segment(self) -> Optional[str]:
        if not self.setup_commands:
            return None

        return "\n".join(self.setup_commands)

    def _sanitized_job_name(self, flow_run: FlowRun) -> str:
        return _JOB_NAME_SANITIZE_PATTERN.sub("_", flow_run.name or str(flow_run.id))


class SlurmCliWorkerTemplateVariables(BaseVariables):
    cpu: int = Field(default=1, description="CPU count required for the flow")
    memory: int = Field(default=4, description="Memory in GB required for the flow")
    partition: Optional[str] = Field(default=None, description="Slurm partition to use")
    account: Optional[str] = Field(default=None, description="Slurm account to charge the job to")
    qos: Optional[str] = Field(default=None, description="Slurm QOS to submit the job under")
    time_limit: int = Field(
        default=1,
        title="Time Limit",
        description="Max wall time in hours for the flow",
    )
    working_dir: Path = Field(
        title="Working Directory",
        description="Directory to run the flow from",
    )
    setup_commands: List[str] = Field(
        default_factory=list,
        title="Setup Commands",
        description=(
            "Shell commands run before the flow-run command, e.g. "
            "'module load nextflow/24.10.3', 'source /path/to/asm_venv/bin/activate'."
        ),
        examples=[["module load nextflow/24.10.3", "source ~/asm_venv/bin/activate"]],
    )
    shebang: str = Field(
        default="#!/bin/bash",
        pattern=r"^#!/.+$",
        description="Indicates which shell to use when running the slurm job",
    )
