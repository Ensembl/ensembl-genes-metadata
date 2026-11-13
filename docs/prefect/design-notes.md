# Design notes

Why `gb_prefect` works the way it does. Each of these replaced something that caused a real
problem; keep them in mind before simplifying.

**Each flow run is one Slurm job (`slurm-cli` worker); tasks don't submit their own jobs.**
The first flows wrapped each Nextflow run in its own `sbatch --wait`, inside a flow running elsewhere:
two layers of Slurm, logs split between them, and resources set in code. Now the worker submits the
flow run itself as the job, the job's resources are deployment settings, and the flow's tasks run
inside it.

**Nextflow and scripts run without a shell.** Tasks run the argument list directly
(`nextflow/runner.py`). Nothing in a parameter can be read as shell syntax, there is no quoting to
get wrong, and no login shell whose profile changes the environment. The environment comes only from
the work pool's `setup_commands`. The `*_command.sh` file is a *record* of the same command, for
reading and re-running by hand.

**Credentials go in a params file or the environment, never the command line.** Command lines are
visible to other users in the process list, end up in Nextflow's `.nextflow.log`, and were being
written to logs and artifacts. A pipeline gets secrets through an owner-only `-params-file`; a script
through `secret_env`. Their values are redacted from everything `gb_prefect` writes.

**A run fails when Nextflow ignored failed tasks.** With `errorStrategy 'ignore'`, a pipeline can
report "Successful" and exit 0 when its key tasks failed; in a BUSCO test, the run directory was
then cleaned up as a success and the GCA left `in_progress` forever. The task now reads Nextflow's
output for `-- Error is ignored` and fails the run (`NextflowRun.fail_on_ignored_errors`).

**Pipelines, scripts and maps are read from ENSCODE.** What runs is the ENSCODE checkout, the same
for the pipeline, its map and its scripts, so they always match each other. Flow code comes from
the deployment's branch, and `gb_metadata` from the venv: three sources, listed in
[Operations](operations.md).

**One spec per pipeline mode, plain Python.** A spec only maps inputs to a `NextflowRun`, with no
Prefect, so it is tested without a server, reused by any flow, and the only place that knows a
pipeline's command line. A pipeline with modes (the statistics pipeline: BUSCO, OMArk, ...) gets a
spec, and a map, per mode, so each run shows only its own steps.

**One deploy script for every deployment.** `deploy_flows.py` describes all of them in code, so a
deployment can be recreated exactly, reviewed in a diff, and changed in one place. Pool names are
defaults overridable per pool, for test pools.

**Live maps never fail a run.** The nf-metro server is a convenience. If it is down, the token is
wrong or the map is missing, the run logs a warning and carries on; with the `nf-metro-url` Variable
unset, nothing is attempted.

**Dry runs need nothing but Python.** No credentials, no ENSCODE, no input files on the machine, no
server: so the command a flow would run can always be checked first.

**Each rerun appends to the same log, and each artifact shows only its own run.** A rerun in the same
folder keeps the history in the log file, while the artifact (what you see in the UI) stays about
the run you're looking at. On a Prefect retry, Nextflow is launched with `-resume`, reusing the
cached tasks in that folder (the flow needs `retries` set on its `@flow` for Prefect to retry it).
