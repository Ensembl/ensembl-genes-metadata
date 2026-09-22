import re
from datetime import datetime
from pathlib import Path
from typing import Optional

from prefect import task  # type: ignore
from prefect.states import Failed, Completed  # type: ignore

from gb_prefect.models.pipeline_options import PipelineCredentials, TaskRunOptions
from gb_prefect.utils.artifact_utils import create_registry_run_artifact
from gb_prefect.utils.enscode_utils import resolve_enscode
from gb_prefect.utils.logging_utils import append_log
from gb_prefect.utils.shell_utils import run_cmd_bash_capture


@task(log_prints=True)
def update_assemblies(
    gca_list: str,
    outdir: str,
    asm_venv: str,
    credentials: PipelineCredentials,
    slack_report: bool = True,
    date: Optional[str] = None,
    enscode: Optional[str] = None,
    run_options: Optional[TaskRunOptions] = None,
):
    """Build and submit the SLURM job that runs the assembly metadata update Nextflow pipeline."""
    if slack_report and not credentials.slack_params:
        raise ValueError("credentials.slack_params is required when slack_report is True")

    run_options = run_options or TaskRunOptions()
    date = date or datetime.now().strftime("%Y-%m-%d")
    outdir_path = Path(outdir)
    log = outdir_path / f"log_flow_update_assemblies_{date}.log"
    command_file = outdir_path / f"update_assemblies_command_{date}.sh"
    log.parent.mkdir(parents=True, exist_ok=True)

    enscode = resolve_enscode(enscode, run_options.dry_run)
    append_log(log, f"[{datetime.now()}] INFO: ENSCODE set to {enscode}.\n")

    slack_params_flag = (
        f"--slack_params '{credentials.slack_params}' \\\n    " if credentials.slack_params else ""
    )

    sbatch_script = f"""#!/bin/bash
#SBATCH --job-name=asm_registry_update_{date}
#SBATCH --output={outdir}/asm_registry_update_slurm_%j.out
#SBATCH --error={outdir}/asm_registry_update_slurm_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=4G
umask 002

module load nextflow/24.10.3
source {asm_venv}/bin/activate

cd {outdir}

nextflow -C {enscode}/ensembl-genes-metadata/pipelines/assembly_metadata_update/nextflow.config \
run {enscode}/ensembl-genes-metadata/pipelines/assembly_metadata_update/main.nf \
    --output_dir {outdir} \
    --gca_list {gca_list} \
    --gca_input true \
    --metadata_params_string '{credentials.metadata_params_string}' \
    --slack_report {str(slack_report).lower()} \
    {slack_params_flag}-with-report \
    -with-dag {outdir}/assembly_update_dag_{date}.png
"""

    append_log(log, f"[{datetime.now()}] INFO: sbatch script:\n{sbatch_script}\n")
    command_file.write_text(sbatch_script)
    command_file.chmod(0o755)

    rc, job_id = _submit_and_collect(command_file, outdir_path, log, run_options.dry_run)

    append_log(log, f"[{datetime.now()}] INFO: Return code {rc}.\n")

    if run_options.create_artifact:
        create_registry_run_artifact(
            date=date,
            outdir=outdir,
            command_file=str(command_file),
            cmd=sbatch_script,
            log_text=log.read_text(),
            rc=rc,
            dry_run=run_options.dry_run,
        )

    result = {
        "returncode": rc,
        "command": sbatch_script,
        "command_file": str(command_file),
        "log_file": str(log),
        "slurm_job_id": job_id if not run_options.dry_run else None,
        "pipeline_run_date": date,
        "pipeline_ran": not run_options.dry_run,
        "dry_run": run_options.dry_run,
    }

    if rc != 0:
        return Failed(
            message=f"Nextflow pipeline failed with return code {rc}",
            data=result,
        )
    return Completed(data=result)


@task(log_prints=True)
def update_assemblies_by_date(
    screen_date: str,
    outdir: str,
    asm_venv: str,
    credentials: PipelineCredentials,
    slack_report: bool = True,
    enscode: Optional[str] = None,
    run_options: Optional[TaskRunOptions] = None,
):
    """Build and submit the SLURM job that screens the metadata database for high-priority
    assemblies released after screen_date and runs the assembly metadata update Nextflow
    pipeline against them (--screen_date mode, mutually exclusive with --gca_input/--gca_list).
    """
    if slack_report and not credentials.slack_params:
        raise ValueError("credentials.slack_params is required when slack_report is True")

    run_options = run_options or TaskRunOptions()
    screen_date = datetime.strptime(screen_date, "%Y-%m-%d").strftime("%Y-%m-%d")
    outdir_path = Path(outdir)
    log = outdir_path / f"log_flow_update_assemblies_by_date_{screen_date}.log"
    command_file = outdir_path / f"update_assemblies_by_date_command_{screen_date}.sh"
    log.parent.mkdir(parents=True, exist_ok=True)

    enscode = resolve_enscode(enscode, run_options.dry_run)
    append_log(log, f"[{datetime.now()}] INFO: ENSCODE set to {enscode}.\n")

    slack_params_flag = (
        f"--slack_params '{credentials.slack_params}' \\\n    " if credentials.slack_params else ""
    )

    sbatch_script = f"""#!/bin/bash
#SBATCH --job-name=asm_registry_update_screen_{screen_date}
#SBATCH --output={outdir}/asm_registry_update_slurm_%j.out
#SBATCH --error={outdir}/asm_registry_update_slurm_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=4G
umask 002

module load nextflow/24.10.3
source {asm_venv}/bin/activate

cd {outdir}

nextflow -C {enscode}/ensembl-genes-metadata/pipelines/assembly_metadata_update/nextflow.config \
run {enscode}/ensembl-genes-metadata/pipelines/assembly_metadata_update/main.nf \
    --output_dir {outdir} \
    --screen_date {screen_date} \
    --metadata_params_string '{credentials.metadata_params_string}' \
    --slack_report {str(slack_report).lower()} \
    {slack_params_flag}-with-report \
    -with-dag {outdir}/assembly_update_dag_{screen_date}.png
"""

    append_log(log, f"[{datetime.now()}] INFO: sbatch script:\n{sbatch_script}\n")
    command_file.write_text(sbatch_script)
    command_file.chmod(0o755)

    rc, job_id = _submit_and_collect(command_file, outdir_path, log, run_options.dry_run)

    append_log(log, f"[{datetime.now()}] INFO: Return code {rc}.\n")

    if run_options.create_artifact:
        create_registry_run_artifact(
            date=screen_date,
            outdir=outdir,
            command_file=str(command_file),
            cmd=sbatch_script,
            log_text=log.read_text(),
            rc=rc,
            dry_run=run_options.dry_run,
        )

    result = {
        "returncode": rc,
        "command": sbatch_script,
        "command_file": str(command_file),
        "log_file": str(log),
        "slurm_job_id": job_id if not run_options.dry_run else None,
        "pipeline_run_date": screen_date,
        "pipeline_ran": not run_options.dry_run,
        "dry_run": run_options.dry_run,
    }

    if rc != 0:
        return Failed(
            message=f"Nextflow pipeline failed with return code {rc}",
            data=result,
        )
    return Completed(data=result)


def _submit_and_collect(command_file: Path, outdir_path: Path, log: Path, dry_run: bool):
    """Submit the sbatch script (unless dry_run) and append its SLURM output/error to the log."""
    if dry_run:
        append_log(log, f"[{datetime.now()}] INFO: Dry run enabled; sbatch job was not submitted.\n")
        return 0, None

    sbatch_result = run_cmd_bash_capture(f"sbatch --wait {command_file}", log_path=log)
    rc = sbatch_result.returncode

    job_id_match = re.search(r"Submitted batch job (\d+)", sbatch_result.stdout)
    job_id = job_id_match.group(1) if job_id_match else None
    append_log(log, f"[{datetime.now()}] INFO: SLURM job ID: {job_id}.\n")

    if job_id:
        for slurm_file in (
            outdir_path / f"asm_registry_update_slurm_{job_id}.out",
            outdir_path / f"asm_registry_update_slurm_{job_id}.err",
        ):
            if slurm_file.exists():
                append_log(log, f"[{datetime.now()}] INFO: --- {slurm_file.name} ---\n")
                append_log(log, slurm_file.read_text())

    return rc, job_id
