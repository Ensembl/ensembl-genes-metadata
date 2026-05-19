from prefect import task # type: ignore
from gb_prefect.utils.logging_utils import append_log 
from gb_prefect.utils.shell_utils import run_cmd_bash_capture
from gb_prefect.utils.artifact_utils import create_registry_run_artifact
from pathlib import Path
from datetime import datetime
import os

@task(log_prints=True)
def register_assemblies(
    date:str,
    outdir: str,
    enscode: str = None,
    dry_run: bool = False,
    create_artifact: bool = True,
):


    outdir_path = Path(outdir)
    date_fmt = datetime.strptime(date, "%m-%d-%Y").strftime("%Y-%m-%d")
    log = outdir_path / f"log_flow_register_assemblies_{date_fmt}.log"
    command_file = outdir_path / f"register_assemblies_command_{date_fmt}.sh"
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
    nextflow run \
    {enscode}/ensembl-genes-metadata/pipelines/assembly_metadata/main.nf \
        --output_dir {outdir} \
        --enscode {enscode} \
        --date {date}
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
        create_registry_run_artifact(
            date=date_fmt,
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
        "log_file": str(log),
        "pipeline_run_date": date_fmt,
        "pipeline_ran": not dry_run,
        "dry_run": dry_run,
    }