import re
from datetime import datetime
from pathlib import Path

from prefect import task  # type: ignore

from gb_prefect.utils.enscode_utils import resolve_enscode
from gb_prefect.utils.logging_utils import append_log
from gb_prefect.utils.shell_utils import run_cmd_bash_capture


@task(log_prints=True)
def is_reference(
    file_path: str,
    output_path: str,
    enscode: str,
    asm_venv: str,
):
    """Build and submit the SLURM job that runs the is_reference check script."""
    outdir_path = Path(output_path)
    output_file = outdir_path / f"{Path(file_path).stem}_output.csv"
    log = outdir_path / f"log_flow_is_reference_{datetime.now().strftime('%Y-%m-%d')}.log"
    command_file = outdir_path / f"is_reference_command_{datetime.now().strftime('%Y-%m-%d')}.sh"
    output_file.parent.mkdir(parents=True, exist_ok=True)
    log.parent.mkdir(parents=True, exist_ok=True)

    enscode = resolve_enscode(enscode, dry_run=False)
    append_log(log, f"[{datetime.now()}] INFO: ENSCODE set to {enscode}.\n")

    sbatch_script = f"""#!/bin/bash
#SBATCH --job-name=is_reference
#SBATCH --output={outdir_path}/slurm_%j.out
#SBATCH --error={outdir_path}/slurm_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=1G

source {asm_venv}/bin/activate
cd {outdir_path}

python {enscode}/ensembl-genes-metadata/src/python/is_reference.py \
    --file {file_path} --output {output_file}
"""
    append_log(log, f"[{datetime.now()}] INFO: sbatch script:\n{sbatch_script}\n")
    command_file.write_text(sbatch_script)
    command_file.chmod(0o755)

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

    if rc != 0:
        raise RuntimeError(f"SLURM job failed with return code {rc}. Check logs for details.")

    append_log(log, f"[{datetime.now()}] INFO: SLURM job completed successfully.\n")
    return str(output_file)
