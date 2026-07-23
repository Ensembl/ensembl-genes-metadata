import re
from datetime import datetime
from pathlib import Path
from typing import Optional

from prefect import task  # type: ignore

from gb_prefect.utils.artifact_utils import create_busco_run_artifact
from gb_prefect.utils.enscode_utils import resolve_enscode
from gb_prefect.utils.logging_utils import append_log
from gb_prefect.utils.shell_utils import run_cmd_bash_capture


@task(log_prints=True)
def run_nextflow_busco(
    csv_file: str,
    outdir: str,
    asm_venv: str,
    enscode: str = None,
    dry_run: bool = False,
    create_artifact: bool = True,
):

    date = datetime.now().strftime("%Y-%m-%d")
    outdir_path = Path(outdir)
    csv_stem = Path(csv_file).stem
    log = outdir_path / f"log_flow_busco_{csv_stem}_{date}.log"
    command_file = outdir_path / f"busco_genome_nextflow_command_{csv_stem}_{date}.sh"
    log.parent.mkdir(parents=True, exist_ok=True)

    enscode = enscode or os.environ.get("ENSCODE")
    if not enscode:
        if dry_run:
            enscode = "<ENSCODE>"
        else:
            raise ValueError(
                "ENSCODE is required when dry_run=False. Pass enscode=..., set ENSCODE, or run with dry_run=True."
            )
    append_log(log, f"[{datetime.now()}] INFO: ENSCODE set to {enscode}.\n")

    sbatch_script = f"""#!/bin/bash
#SBATCH --job-name=busco_genome_{csv_stem}_{date}
#SBATCH --output={outdir}/busco_genome_slurm_%j.out
#SBATCH --error={outdir}/busco_genome_slurm_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=2G
umask 002

module load nextflow/24.10.3
source {asm_venv}/bin/activate

cd {outdir}
nextflow run {enscode}/ensembl-genes-nf/pipelines/statistics/main.nf \
    --csvFile {csv_file} \
    --run_busco_ncbi \
    --outdir {outdir}
"""
    append_log(log, f"[{datetime.now()}] INFO: sbatch script:\n{sbatch_script}\n")
    command_file.write_text(sbatch_script)
    command_file.chmod(0o755)

    if dry_run:
        rc = 0
        job_id = None
        append_log(log, f"[{datetime.now()}] INFO: Dry run enabled; Nextflow command was not executed.\n")
    else:
        sbatch_result = run_cmd_bash_capture(f"sbatch --wait {command_file}", log_path=log)
        rc = sbatch_result.returncode

        job_id_match = re.search(r"Submitted batch job (\d+)", sbatch_result.stdout)
        job_id = job_id_match.group(1) if job_id_match else None
        append_log(log, f"[{datetime.now()}] INFO: SLURM job ID: {job_id}.\n")

        if job_id:
            for slurm_file in (
                outdir_path / f"busco_genome_slurm_{job_id}.out",
                outdir_path / f"busco_genome_slurm_{job_id}.err",
            ):
                if slurm_file.exists():
                    append_log(log, f"[{datetime.now()}] INFO: --- {slurm_file.name} ---\n")
                    append_log(log, slurm_file.read_text())

    append_log(log, f"[{datetime.now()}] INFO: Return code {rc}.\n")

    if create_artifact:
        create_busco_run_artifact(
            csv_file=csv_file,
            outdir=outdir,
            command_file=str(command_file),
            cmd=sbatch_script,
            log_text=log.read_text(),
            rc=rc,
            dry_run=dry_run,
        )

    return {
        "returncode": rc,
        "command": sbatch_script,
        "command_file": str(command_file),
        "log": str(log),
        "pipeline_ran": not dry_run,
        "dry_run": dry_run,
    }
