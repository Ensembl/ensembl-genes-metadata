from prefect import task  # type: ignore
from gb_prefect.utils.artifact_utils import create_busco_run_artifact
from gb_prefect.utils.logging_utils import append_log
from gb_prefect.utils.shell_utils import run_cmd_bash_capture
from datetime import datetime
import os
from pathlib import Path


@task(log_prints=True)
def run_nextflow_busco(
    csv_file: str,
    outdir: str,
    enscode: str = None,
    dry_run: bool = False,
    create_artifact: bool = True,
):

    outdir_path = Path(outdir)
    csv_stem = Path(csv_file).stem
    log = outdir_path / f"log_flow_busco_{csv_stem}.log"
    command_file = outdir_path / f"busco_genome_nextflow_command_{csv_stem}.sh"
    log.parent.mkdir(parents=True, exist_ok=True)

    enscode = enscode or os.environ.get("ENSCODE")
    if not enscode:
        if dry_run:
            enscode = "<ENSCODE>"
        else:
            raise ValueError("ENSCODE is required when dry_run=False. Pass enscode=..., set ENSCODE, or run with dry_run=True.")
    append_log(log, f"[{datetime.now()}] INFO: ENSCODE set to {enscode}.\n")
    
    cmd = f"""
    module load ensembl/asm_update_dev &&
    module load nextflow &&
    cd {outdir} &&
    nextflow run {enscode}/ensembl-genes-nf/pipelines/statistics/main.nf \
        --csvFile {csv_file} \
        --run_busco_ncbi \
        --outdir {outdir}/busco_genome_output
    """
    append_log(log, f"[{datetime.now()}] INFO: {cmd}.\n")
    command_file.write_text(cmd.strip() + "\n")

    if dry_run:
        rc = 0
        append_log(log, f"[{datetime.now()}] INFO: Dry run enabled; Nextflow command was not executed.\n")
    else:
        result = run_cmd_bash_capture(cmd, log_path=log)
        rc = result.returncode

    append_log(log, f"[{datetime.now()}] INFO: Return code {rc}.\n")

    if create_artifact:
        create_busco_run_artifact(
            csv_file=csv_file,
            outdir=outdir,
            command_file=str(command_file),
            cmd=cmd,
            log_text=log.read_text(),
            rc=rc,
            dry_run=dry_run,
        )
            
    return {
        "returncode": rc,
        "command": cmd,
        "command_file": str(command_file),
        "log": str(log),
        "pipeline_ran": not dry_run,
        "dry_run": dry_run,
    }
