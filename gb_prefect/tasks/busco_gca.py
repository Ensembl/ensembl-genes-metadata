from datetime import datetime
from pathlib import Path
from typing import Optional

from prefect import task  # type: ignore
from prefect.states import Completed  # type: ignore

from gb_prefect.utils.artifact_utils import create_busco_run_artifact
from gb_prefect.utils.enscode_utils import resolve_enscode
from gb_prefect.utils.logging_utils import append_log
from gb_prefect.utils.shell_utils import run_cmd_bash_capture
from gb_prefect.utils.stats_split_csv import write_single_gca_csv


@task(log_prints=True)
def run_nextflow_busco_gca(
    gca: str,
    taxon_id: str,
    outdir: str,
    busco_dataset: Optional[str] = None,
    enscode: Optional[str] = None,
    dry_run: bool = False,
    create_artifact: bool = True,
):
    """Run the BUSCO genome-statistics Nextflow pipeline for a single GCA.

    Unlike gb_prefect.tasks.busco.run_nextflow_busco, this does not build or submit
    its own Slurm job -- it assumes the flow run it belongs to is already executing
    inside a Slurm job submitted by the slurm-cli worker (see gb_prefect/worker/),
    with `nextflow` already on PATH via the worker's `setup_commands`. It just runs
    `nextflow run` directly.
    """
    date = datetime.now().strftime("%Y-%m-%d")
    outdir_path = Path(outdir)
    log = outdir_path / f"log_flow_busco_{gca}_{date}.log"
    command_file = outdir_path / f"busco_gca_nextflow_command_{gca}_{date}.sh"
    log.parent.mkdir(parents=True, exist_ok=True)

    enscode = resolve_enscode(enscode, dry_run)
    append_log(log, f"[{datetime.now()}] INFO: ENSCODE set to {enscode}.\n")

    csv_file = write_single_gca_csv(gca, taxon_id, outdir, busco_dataset)
    append_log(log, f"[{datetime.now()}] INFO: Wrote single-row CSV to {csv_file}.\n")

    dataset_flag = f"--buscoDataset {busco_dataset}" if busco_dataset else ""

    nextflow_command = f"""#!/bin/bash
cd {outdir}
export ENSCODE={enscode}
nextflow run {enscode}/ensembl-genes-nf/pipelines/statistics/main.nf \
    --csvFile {csv_file} \
    --run_busco_ncbi \
    --outdir {outdir} \
    --enscode {enscode} \
    {dataset_flag}
"""
    append_log(log, f"[{datetime.now()}] INFO: Nextflow command:\n{nextflow_command}\n")
    command_file.write_text(nextflow_command)
    command_file.chmod(0o755)

    if dry_run:
        rc = 0
        append_log(log, f"[{datetime.now()}] INFO: Dry run enabled; Nextflow command was not executed.\n")
    else:
        # NOTE: run_cmd_bash_capture() overwrites (not appends to) log_path with the
        # subprocess's own stdout/stderr -- the ENSCODE/CSV/command lines appended above
        # are lost once this runs. Pre-existing behavior shared with gb_prefect.tasks.busco
        # and gb_prefect.tasks.registry_update; tracked separately, not fixed here.
        result = run_cmd_bash_capture(f"bash {command_file}", log_path=log)
        rc = result.returncode

    append_log(log, f"[{datetime.now()}] INFO: Return code {rc}.\n")

    if create_artifact:
        create_busco_run_artifact(
            csv_file=csv_file,
            outdir=outdir,
            command_file=str(command_file),
            cmd=nextflow_command,
            log_text=log.read_text(),
            rc=rc,
            dry_run=dry_run,
        )

    result_data = {
        "returncode": rc,
        "command": nextflow_command,
        "command_file": str(command_file),
        "log": str(log),
        "csv_file": csv_file,
        "gca": gca,
        "pipeline_ran": not dry_run,
        "dry_run": dry_run,
    }

    if rc != 0:
        # Raise rather than `return Failed(data=result_data)`: busco_gca_flow just returns
        # this task's result directly, and a flow propagating a Failed state whose data is a
        # plain dict (not an exception) crashes the engine with "dict cannot be resolved into
        # an exception" during result resolution -- confirmed independent of anything BUSCO-
        # specific. Raising a real exception avoids that; result_data is still fully captured
        # in the log file and the artifact above.
        raise RuntimeError(
            f"Nextflow BUSCO pipeline failed for {gca} with return code {rc}. "
            f"See log: {result_data['log']}"
        )
    return Completed(data=result_data)
