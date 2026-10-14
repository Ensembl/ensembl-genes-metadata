"""
Collection of Prefect Nextflow input and output classes. 
"""

from dataclasses import dataclass, field
from typing import Any, Dict, Mapping, Optional, Sequence


@dataclass(frozen=True)
class CsvInput:
    """Rows the task writes to <outdir>/<run name>.csv before the run.

    Attributes:
        param: Pipeline parameter the CSV path is passed as, e.g. "csvFile".
        columns: Header row, in order.
        rows: One dict per data row, keyed by column; missing keys are written empty.
    """

    param: str
    columns: Sequence[str]
    rows: Sequence[Mapping[str, Any]]


@dataclass(frozen=True)
class NextflowRun:  # pylint: disable=too-many-instance-attributes
    """Everything needed to launch one Nextflow pipeline run.

    Attributes:
        name: Short identifier
        enscode: Resolved ENSCODE directory.
        pipeline: Path of the pipeline's main.nf, relative to enscode.
        outdir: Launch directory of the run
        params: Pipeline parameters passed on the command line.
        secret_params: Sensitive parameters.
        config: Config file passed with -C, relative to enscode. -C makes Nextflow use ONLY this
            file: the pipeline's own nextflow.config is not loaded unless this is that file.
        csv_input: Optional rows written to a CSV and passed as csv_input.param.
        profile: Nextflow profile(s) (-profile).
        with_report: Add -with-report.
        dag_file: Write the DAG to <outdir>/<dag_file> (-with-dag).
        weblog_url: Send run events to nf-metro serve endpoint (-with-weblog).
        resume: Add -resume.
    """

    name: str
    enscode: str
    pipeline: str
    outdir: str
    params: Mapping[str, Any] = field(default_factory=dict)
    secret_params: Mapping[str, Any] = field(default_factory=dict)
    config: Optional[str] = None
    csv_input: Optional[CsvInput] = None
    profile: Optional[str] = None
    with_report: bool = False
    dag_file: Optional[str] = None
    weblog_url: Optional[str] = None
    resume: bool = False

    def __post_init__(self) -> None:
        """Reject a parameter given twice: Nextflow would silently pick one.

        Command-line params override the -params-file, so a key in both params and
        secret_params would quietly drop the secret value; the CSV parameter likewise.
        """
        overlap = set(self.params) & set(self.secret_params)
        if overlap:
            raise ValueError(f"Parameter(s) in both params and secret_params: {sorted(overlap)}")
        if self.csv_input is not None and (
            self.csv_input.param in self.params or self.csv_input.param in self.secret_params
        ):
            raise ValueError(
                f"csv_input.param {self.csv_input.param!r} is also given in params/secret_params"
            )


@dataclass(frozen=True)
class NextflowResult:  # pylint: disable=too-many-instance-attributes
    """Outcome of run_nextflow_pipeline, returned to the calling flow.

    Attributes:
        name: The run's name.
        returncode: Nextflow's exit code (0 for a dry run).
        command: The command script.
        command_file: Path of the command script.
        log_file: Path of the flow log holding the command and Nextflow's output.
        outdir: The run's launch/output directory.
        dry_run: True if Nextflow was not executed.
        resumed: True if the run was launched with -resume.
        files: Other files written for the run (csv, params, ...), by type.
    """

    name: str
    returncode: int
    command: str
    command_file: str
    log_file: str
    outdir: str
    dry_run: bool
    resumed: bool
    files: Dict[str, str] = field(default_factory=dict)


class NextflowRunError(RuntimeError):
    """Raised when a Nextflow run exits non-zero; result holds everything about the run."""

    def __init__(self, result: NextflowResult):
        super().__init__(
            f"Nextflow run {result.name} failed with return code {result.returncode}. "
            f"See log: {result.log_file}"
        )
        self.result = result

    def __reduce__(self):
        # Default exception pickling would call cls(message); rebuild from result instead.
        return (self.__class__, (self.result,))
