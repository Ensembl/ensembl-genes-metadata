import shutil
from pathlib import Path
from typing import Optional

from prefect import task  # type: ignore

from gb_metadata.update_busco_events import mark_busco_failed
from gb_prefect.utils.credentials_utils import PipelineCredentials, resolve_credentials

# Post-run steps for flows running ensembl-genes-nf's statistics pipeline, one run directory per
# key (e.g. <outdir>/<gca>): cleaning that directory is shared by every statistics mode;
# record_busco_failure is BUSCO's own status update.

FAILED_RUN_CLEANUP = ["work", "cache", ".nextflow", "nextflow_*_params.json"]


@task(log_prints=True)
def record_busco_failure(
    gca: str,
    credentials: Optional[PipelineCredentials] = None,
    metadata_secret_block: Optional[str] = None,
) -> None:
    """Mark genome_busco.status as failed for gca after its BUSCO run failed."""
    credentials = resolve_credentials(credentials, metadata_secret_block=metadata_secret_block)
    mark_busco_failed(gca, credentials.metadata_db_params, execute=True)
    print(f"Marked {gca} as failed")


@task(log_prints=True)
def cleanup_statistics_outdir(run_outdir: str, success: bool) -> None:
    """Clean a statistics run directory (<outdir>/<key>) after an automatic run.

    success=True removes the whole directory.
    success=False removes only what matches FAILED_RUN_CLEANUP, keeping .nextflow.log,
    the flow log and the command file.
    """
    run_path = Path(run_outdir)
    if success:
        targets = [run_path]
    else:
        targets = sorted({target for pattern in FAILED_RUN_CLEANUP for target in run_path.glob(pattern)})

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
