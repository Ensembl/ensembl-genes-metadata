import subprocess
from datetime import datetime
from pathlib import Path

from prefect import get_run_logger, task # type: ignore
from prefect.artifacts import create_markdown_artifact  # type: ignore

from gb_prefect.models.asm_outcomes import AssemblyOutcome
from gb_prefect.models.paths import RunPaths
from gb_prefect.utils.logging_utils import append_log, filter_stderr
from gb_prefect.utils.shell_utils import run_cmd_bash_capture
from gb_prefect.utils.logging_utils import append_log



@task
def run_assembly_pipeline(
    paths: RunPaths, 
    update_date: str,
    enscode: str = None,
    dry_run: bool = False,
    create_artifact: bool = True,
) -> AssemblyOutcome:
    
    # Initial set up
    logger = get_run_logger()
    run_dir = Path(paths["run_dir"])
    cron_log = Path(paths["cron_log"])
    asm_log = Path(paths["asm_log"])
    out_dir = run_dir / "nextflow output"

    if dry_run:
        # Create a placeholder for ENSCODE
        enscode = "<ENSCODE>"
    else:
        # Only attempt to load modules when not in dry run mode
        check_env = (
            "module load ensembl/asm_update_dev && "
            "module load nextflow && "
            'if [[ -z "${ENSCODE:-}" ]]; then echo "ENSCODE_NOT_SET"; else echo "$ENSCODE"; fi'
        )
        env_proc = subprocess.run(
            ["bash", "-lc", f"source /etc/bashrc && {check_env}"],
            text=True,
            capture_output=True,
        )
        stderr_clean = filter_stderr(env_proc.stderr or "")
        stdout_lines = [l.strip() for l in (env_proc.stdout or "").splitlines() if l.strip()]
        enscode = stdout_lines[-1] if stdout_lines else ""

        if env_proc.returncode != 0 or enscode == "ENSCODE_NOT_SET" or not enscode:
            append_log(cron_log, f"[{datetime.now()}] ERROR: ENSCODE not set after module load.\n")
            if stderr_clean:
                append_log(cron_log, f"[{datetime.now()}] stderr (filtered):\n{stderr_clean}\n")

        return {
            "option": "option4",
            "returncode": env_proc.returncode if env_proc.returncode != 0 else 99,
            "run_dir": str(run_dir),
            "out_dir": str(out_dir),
            "assemblies_to_register": str(out_dir / "assemblies_to_register.txt"),
            "gca_to_run_ncbi": str(out_dir / "gca_to_run_ncbi.csv"),
            "report": str(out_dir / "report.txt"),
            "assemblies_nonempty": False,
            "gca_exists_nonempty": False,
        }

    append_log(cron_log, f"[{datetime.now()}] ENSCODE is set to: {enscode}\n")
    logger.info(f"ENSCODE={enscode}")

    assembly_cmd = (
        "source /var/tmp/prefect_venv/asm_venv/bin/activate && "
        f"cd {run_dir.parent} && "
        "nextflow run "
        f"{enscode}/ensembl-genes-metadata/pipelines/assembly_metadata/main.nf "
        f"--output_dir {run_dir} "
        f"--enscode {enscode} "
        f"--date {update_date} "
    )

    append_log(cron_log, f"[{datetime.now()}] Running Assembly Pipeline:\n{assembly_cmd}\n")

    if dry_run:
        rc=0
        append_log(cron_log, f"[{datetime.now()}] Dry run enabled; Assembly pipeline command was not executed.\n")
    else:
        result = run_cmd_bash_capture(assembly_cmd, log_path=asm_log)
        rc = result.returncode

    
    assemblies = out_dir / "assemblies_to_register.txt"
    gca_csv = out_dir / "gca_to_run_ncbi.csv"
    report = out_dir / "report.txt"

    assemblies.parent.mkdir(parents=True, exist_ok=True)
    if not assemblies.exists():
        assemblies.write_text("")

    assemblies_nonempty = assemblies.stat().st_size > 0
    gca_exists_nonempty = gca_csv.exists() and gca_csv.stat().st_size > 0

    if rc == 0:
        if gca_exists_nonempty and report.exists():
            option = "option1" # run succedded with assemblies to run BUSCO
        else:
            option = "option2" # run succeeded but no assemblies to run  BUSCO"
    else:
        if assemblies_nonempty:
            option = "option3" # Run started but failed before completion
        else:
            option = "option4" # Run failded but to environment or module issues

    append_log(
        cron_log,
        f"[{datetime.now()}] Assembly pipeline returncode={rc}, classified as {option}\n",
    )

    if result.stderr:
        filtered = filter_stderr(result.stderr)
        if filtered:
            append_log(cron_log, f"[{datetime.now()}] stderr (filtered):\n{filtered}\n")

        if create_artifact:
        log_text = log.read_text()
        create_markdown_artifact(
            key="busco-run-log",
            description="Nextflow BUSCO run log",
            markdown=f"""# BUSCO run log

            **CSV file:** `{csv_file}`  
            **Outdir:** `{outdir}`  
            **Dry run:** `{dry_run}`
            **Command file:** `{command_file}`
            **Return code:** `{rc}`

            ## Nextflow command

            ```bash
            {cmd}
            ```

            ## Log output

            ```text
            {log_text}
            ```
        """
        )

    return {
        "option": option,
        "returncode": rc,
        "run_dir": str(run_dir),
        "out_dir": str(out_dir),
        "assemblies_to_register": str(assemblies),
        "gca_to_run_ncbi": str(gca_csv),
        "report": str(report),
        "assemblies_nonempty": assemblies_nonempty,
        "gca_exists_nonempty": gca_exists_nonempty,
    }