"""Input and output of run_script (gb_prefect/tasks/script.py): one command, e.g. a Python script.

The script counterpart of models/nextflow.py: a flow builds a ScriptRun, the task runs it inside
the flow run's own Slurm job, and returns a ScriptResult or raises ScriptRunError.
"""

from dataclasses import dataclass, field
from typing import Mapping, Sequence


@dataclass(frozen=True)
class ScriptRun:
    """One command to run, e.g. ``[sys.executable, ".../src/python/is_reference.py", "--file", path]``.

    Attributes:
        name: Short identifier; names the log, command record and artifact, so it must be unique
            within outdir.
        argv: The command and its arguments, run directly (no shell): no quoting needed, and no
            shell syntax (pipes, globs, $VARS) is interpreted.
        outdir: Directory the command runs from; the log and command record are written there.
        env: Extra environment variables, on top of the flow run's own environment (the work
            pool's setup_commands: modules, venv).
        secret_env: Extra environment variables holding secrets (passwords, tokens): their values
            are redacted from the log and artifact and never written to the command record.
            Pass credentials this way, not in argv.
    """

    name: str
    argv: Sequence[str]
    outdir: str
    env: Mapping[str, str] = field(default_factory=dict)
    secret_env: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.argv:
            raise ValueError("argv must name a command to run")
        overlap = set(self.env) & set(self.secret_env)
        if overlap:
            raise ValueError(f"Variable(s) in both env and secret_env: {sorted(overlap)}")


@dataclass(frozen=True)
class ScriptResult:
    """Outcome of run_script, returned to the calling flow.

    Attributes:
        name: The run's name.
        returncode: The command's exit code (0 for a dry run).
        command: The command record (secret_env values left out).
        command_file: Path of the command record.
        log_file: Path of the log holding the command and its output.
        outdir: The run's directory.
        dry_run: True if the command was not executed.
    """

    name: str
    returncode: int
    command: str
    command_file: str
    log_file: str
    outdir: str
    dry_run: bool


class ScriptRunError(RuntimeError):
    """Raised when a script exits non-zero; result holds everything about the run."""

    def __init__(self, result: ScriptResult):
        super().__init__(
            f"Script run {result.name} failed with return code {result.returncode}. "
            f"See log: {result.log_file}"
        )
        self.result = result

    def __reduce__(self):
        # Default exception pickling would call cls(message); rebuild from result instead.
        return (self.__class__, (self.result,))
