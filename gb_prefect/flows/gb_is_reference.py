
from prefect import flow  # type: ignore
from gb_prefect.tasks.is_reference import is_reference


@flow(name="gb_is_reference", log_prints=True)
def gb_is_reference_flow(
    file_path: str,
    output_path: str,
    enscode: str,
    asm_venv: str,
):
    """Run the is_reference check script for file_path via a SLURM job."""
    return is_reference(
        file_path=file_path,
        output_path=output_path,
        enscode=enscode,
        asm_venv=asm_venv,
    )
