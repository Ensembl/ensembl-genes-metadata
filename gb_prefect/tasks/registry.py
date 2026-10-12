import re
from datetime import datetime
from pathlib import Path
from typing import Optional

from prefect import task  # type: ignore
from prefect.states import Failed, Completed  # type: ignore

from gb_prefect.models.pipeline_options import PipelineCredentials, TaskRunOptions
from gb_prefect.utils.artifact_utils import create_registry_run_artifact
from gb_prefect.utils.credentials_utils import redact_credentials
from gb_prefect.utils.enscode_utils import resolve_enscode
from gb_prefect.utils.logging_utils import append_log
from gb_prefect.utils.shell_utils import run_cmd_bash_capture


def _registry_mode_flags(date: Optional[str], gca_list: Optional[str]) -> str:
    """Nextflow flags selecting the assembly_metadata input mode.

    - gca_list: register exactly those accessions (--add_gca); no NCBI screening.
    - date (YYYY-MM-DD): screen NCBI for assemblies released after that date.
    - neither: screen NCBI from the DB's last regular update date.

    The pipeline's --full_screen mode is a developer option and is deliberately not exposed here.
    """
    flags = []
    if gca_list:
        flags.append(f"--add_gca true --gca_list {gca_list}")
    if date:
        flags.append(f"--date {datetime.strptime(date, '%Y-%m-%d').strftime('%m/%d/%Y')}")
    return "".join(f"    {flag} \\\n" for flag in flags)


@task(log_prints=True)
def register_assemblies(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    outdir: str,
    asm_venv: str,
    credentials: PipelineCredentials,
    date: Optional[str] = None,
    gca_list: Optional[str] = None,
    enscode: Optional[str] = None,
    run_options: Optional[TaskRunOptions] = None,
):
    """Build and submit the SLURM job that runs the assembly registry Nextflow pipeline."""
    run_options = run_options or TaskRunOptions()
    mode_flags = _registry_mode_flags(date, gca_list)
    run_date = datetime.now().strftime("%Y-%m-%d")
    outdir_path = Path(outdir)
    log = outdir_path / f"log_flow_register_assemblies_{run_date}.log"
    command_file = outdir_path / f"register_assemblies_command_{run_date}.sh"
    log.parent.mkdir(parents=True, exist_ok=True)

    enscode = resolve_enscode(enscode, run_options.dry_run)
    append_log(log, f"[{datetime.now()}] INFO: ENSCODE set to {enscode}.\n")

    sbatch_script = f"""#!/bin/bash
#SBATCH --job-name=asm_registry_{run_date}
#SBATCH --output={outdir}/slurm_%j.out
#SBATCH --error={outdir}/slurm_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=4G
umask 002

module load nextflow/24.10.3
source {asm_venv}/bin/activate

cd {outdir}

nextflow run {enscode}/ensembl-genes-metadata/pipelines/assembly_metadata/main.nf \\
    --output_dir {outdir} \\
    --enscode {enscode} \\
{mode_flags}    --metadata_params_string '{credentials.metadata_params_string}' \\
    -with-report \\
    -with-dag
"""

    # The command file holds the real credentials so it is owner-only; everything that is
    # logged, returned or published as an artifact gets the redacted copy.
    redacted_script = redact_credentials(sbatch_script, credentials)
    append_log(log, f"[{datetime.now()}] INFO: sbatch script:\n{redacted_script}\n")
    command_file.write_text(sbatch_script)
    command_file.chmod(0o700)

    if run_options.dry_run:
        rc = 0
        job_id = None
        append_log(log, f"[{datetime.now()}] INFO: Dry run enabled; sbatch job was not submitted.\n")
    else:
        sbatch_result = run_cmd_bash_capture(f"sbatch --wait {command_file}", log_path=log)
        rc = sbatch_result.returncode

        job_id_match = re.search(r"Submitted batch job (\d+)", sbatch_result.stdout)
        job_id = job_id_match.group(1) if job_id_match else None
        append_log(log, f"[{datetime.now()}] INFO: SLURM job ID: {job_id}.\n")

        if job_id:
            for slurm_file in (
                outdir_path / f"slurm_{job_id}.out",
                outdir_path / f"slurm_{job_id}.err",
            ):
                if slurm_file.exists():
                    append_log(log, f"[{datetime.now()}] INFO: --- {slurm_file.name} ---\n")
                    append_log(log, slurm_file.read_text())

    append_log(log, f"[{datetime.now()}] INFO: Return code {rc}.\n")

    if run_options.create_artifact:
        create_registry_run_artifact(
            date=run_date,
            outdir=outdir,
            command_file=str(command_file),
            cmd=redacted_script,
            log_text=redact_credentials(log.read_text(), credentials),
            rc=rc,
            dry_run=run_options.dry_run,
        )

    result = {
        "returncode": rc,
        "command": redacted_script,
        "command_file": str(command_file),
        "log_file": str(log),
        "slurm_job_id": job_id if not run_options.dry_run else None,
        "pipeline_run_date": run_date,
        "pipeline_ran": not run_options.dry_run,
        "dry_run": run_options.dry_run,
    }

    if rc != 0:
        return Failed(
            message=f"Nextflow pipeline failed with return code {rc}",
            data=result,
        )
    return Completed(data=result)
