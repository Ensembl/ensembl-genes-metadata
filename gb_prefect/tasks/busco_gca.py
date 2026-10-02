import json
import shutil
from datetime import datetime
from pathlib import Path
from typing import Optional

from prefect import task  # type: ignore
from prefect.states import Completed  # type: ignore

from gb_metadata.update_busco_events import mark_busco_failed
from gb_prefect.models.pipeline_options import PipelineCredentials
from gb_prefect.utils.artifact_utils import create_busco_run_artifact
from gb_prefect.utils.credentials_utils import (
    DEFAULT_METADATA_SECRET_BLOCK,
    DEFAULT_SLACK_SECRET_BLOCK,
    resolve_credentials,
)
from gb_prefect.utils.enscode_utils import resolve_enscode
from gb_prefect.utils.logging_utils import append_log
from gb_prefect.utils.shell_utils import run_cmd_bash_capture
from gb_prefect.utils.stats_split_csv import write_single_gca_csv


@task(log_prints=True)
def run_nextflow_busco_gca(  # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals
    gca: str,
    taxon_id: str,
    outdir: str,
    busco_dataset: Optional[str] = None,
    enscode: Optional[str] = None,
    dry_run: bool = False,
    create_artifact: bool = True,
    update_registry: bool = False,
    credentials: Optional[PipelineCredentials] = None,
    metadata_secret_block: str = DEFAULT_METADATA_SECRET_BLOCK,
):
    """Run the BUSCO genome-statistics Nextflow pipeline for a single GCA.

    Unlike gb_prefect.tasks.busco.run_nextflow_busco, this does not build or submit
    its own Slurm job -- it assumes the flow run it belongs to is already executing
    inside a Slurm job submitted by the slurm-cli worker (see gb_prefect/worker/),
    with `nextflow` already on PATH via the worker's `setup_commands`. It just runs
    `nextflow run` directly.

    update_registry=True adds the pipeline's --update_registry, which loads the genome
    results into the assembly metadata DB and marks genome_busco.status as done. The DB
    connection params it needs (--asm_metadata) are written to a params file readable only
    by the owner and passed with -params-file, so they never appear in the command file, the
    log or the artifact. credentials is optional: when omitted (the normal case for a
    deployment trigger), they're loaded from the metadata_secret_block Prefect Secret block.
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

    registry_flag = ""
    if update_registry:
        params_file = outdir_path / "asm_metadata_params.json"
        registry_flag = f"--update_registry -params-file {params_file}"
        if not dry_run:
            credentials = resolve_credentials(
                credentials, metadata_secret_block, DEFAULT_SLACK_SECRET_BLOCK, slack_report=False
            )
            # Create with owner-only permissions before writing, so the DB password is never
            # readable by anyone else, even briefly.
            params_file.touch(mode=0o600, exist_ok=True)
            params_file.chmod(0o600)
            params_file.write_text(
                json.dumps({"asm_metadata": json.loads(credentials.metadata_params_string)})
            )
            append_log(log, f"[{datetime.now()}] INFO: Wrote asm_metadata params file {params_file}.\n")

    nextflow_command = f"""#!/bin/bash
cd {outdir}
export ENSCODE={enscode}
nextflow run {enscode}/ensembl-genes-nf/pipelines/statistics/main.nf \
    --csvFile {csv_file} \
    --run_busco_ncbi \
    --outdir {outdir} \
    --enscode {enscode} \
    {dataset_flag} \
    {registry_flag}
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
        "update_registry": update_registry,
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


# Nextflow scratch plus the DB credentials file -- removed after a failed automatic run,
# while the logs and command script next to them are kept for debugging.
FAILED_RUN_CLEANUP = ["work", "cache", ".nextflow", "asm_metadata_params.json"]


@task(log_prints=True)
def record_busco_failure(
    gca: str,
    credentials: Optional[PipelineCredentials] = None,
    metadata_secret_block: str = DEFAULT_METADATA_SECRET_BLOCK,
) -> None:
    """Mark genome_busco.status as failed for gca after its BUSCO run failed."""
    credentials = resolve_credentials(
        credentials, metadata_secret_block, DEFAULT_SLACK_SECRET_BLOCK, slack_report=False
    )
    mark_busco_failed(gca, json.loads(credentials.metadata_params_string), execute=True)
    print(f"Marked {gca} as failed")


@task(log_prints=True)
def cleanup_busco_outdir(run_outdir: str, success: bool) -> None:
    """Clean a BUSCO_gca run directory (<outdir>/<gca>) after an automatic run.

    success=True removes the whole directory -- the results are already loaded into the
    assembly metadata DB by the pipeline's --update_registry. success=False removes only
    FAILED_RUN_CLEANUP, keeping .nextflow.log, the flow log and the command script.

    Removal errors are reported but not raised: a cleanup problem shouldn't turn a
    successful BUSCO run into a failed one (or hide the original error of a failed one).
    """
    run_path = Path(run_outdir)
    targets = [run_path] if success else [run_path / name for name in FAILED_RUN_CLEANUP]

    for target in targets:
        if not target.exists():
            continue
        try:
            if target.is_dir():
                shutil.rmtree(target)
            else:
                target.unlink()
            print(f"Removed {target}")
        except OSError as err:
            print(f"WARNING: could not remove {target}: {err}")
