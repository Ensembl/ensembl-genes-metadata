"""
Collection of Prefect Nextflow input and output classes. 
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Mapping, Optional, Sequence


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
        metro_map: The pipeline's nf-metro map (.mmd), relative to enscode. When an nf-metro
            server is configured, the task registers the run with it and sets weblog_url, for a
            live map of the run (see gb_prefect/nextflow/metro.py).
        resume: Add -resume.
        fail_on_ignored_errors: Fail the run when Nextflow exits 0 but its output shows tasks that
            failed and were ignored (errorStrategy 'ignore').
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
    metro_map: Optional[str] = None
    resume: bool = False
    fail_on_ignored_errors: bool = True

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
        ignored_errors: Tasks that failed and were ignored by the pipeline, as
            "<task> (exit <status>)".
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
    ignored_errors: List[str] = field(default_factory=list)


class NextflowRunError(RuntimeError):
    """Raised when a Nextflow run exits non-zero, or exits 0 with ignored task failures (see
    NextflowRun.fail_on_ignored_errors); result holds everything about the run."""

    def __init__(self, result: NextflowResult):
        if result.returncode != 0:
            reason = f"failed with return code {result.returncode}"
        else:
            shown = ", ".join(result.ignored_errors[:5])
            more = f" and {len(result.ignored_errors) - 5} more" if len(result.ignored_errors) > 5 else ""
            reason = (
                f"exited 0, but {len(result.ignored_errors)} task(s) failed and were ignored by the "
                f"pipeline: {shown}{more}"
            )
        super().__init__(f"Nextflow run {result.name} {reason}. See log: {result.log_file}")
        self.result = result

    def __reduce__(self):
        # Default exception pickling would call cls(message); rebuild from result instead.
        return (self.__class__, (self.result,))
