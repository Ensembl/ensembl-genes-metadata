import sys
from datetime import datetime
from pathlib import Path
from typing import Optional

from prefect import flow  # type: ignore

from gb_prefect.models.pipeline_options import TaskRunOptions
from gb_prefect.models.script import ScriptRun
from gb_prefect.tasks.script import run_script
from gb_prefect.utils.enscode_utils import resolve_enscode
from gb_prefect.utils.path_utils import check_input_file, run_folder

IS_REFERENCE_SCRIPT = "ensembl-genes-metadata/src/python/is_reference.py"


@flow(name="gb_is_reference", log_prints=True)
def gb_is_reference_flow(
    file_path: str,
    outdir: str,
    folder_name: Optional[str] = None,
    enscode: Optional[str] = None,
    dry_run: bool = False,
) -> str:
    """Check whether the GCA accessions in file_path are the NCBI reference assembly for their
    taxon, with src/python/is_reference.py; return the path of the CSV it writes.

    file_path: absolute path to a text/CSV file with one GCA accession per line (first column).
    outdir, folder_name: the run goes to <outdir>/<folder_name>; folder_name defaults to
        is_reference_<today>. The CSV is <file_path stem>_output.csv there.
    enscode: ENSCODE directory holding ensembl-genes-metadata; defaults to $ENSCODE.
    dry_run: write and log the command without running it.

    The script runs with run_script, inside this flow run's Slurm job, using the same Python
    (the venv the work pool activates). It needs no credentials: it queries NCBI's public API.
    """
    run_outdir = run_folder(outdir, folder_name, f"is_reference_{datetime.now().strftime('%Y-%m-%d')}")
    output_file = Path(run_outdir) / f"{Path(file_path).stem}_output.csv"

    run = ScriptRun(
        name="is_reference",
        argv=[
            sys.executable,
            f"{resolve_enscode(enscode, dry_run)}/{IS_REFERENCE_SCRIPT}",
            "--file",
            check_input_file(file_path, dry_run),
            "--output",
            str(output_file),
        ],
        outdir=run_outdir,
    )
    run_script(run, TaskRunOptions(dry_run=dry_run))
    return str(output_file)
