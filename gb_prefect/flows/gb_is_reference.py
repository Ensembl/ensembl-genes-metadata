import argparse
import sys
from pathlib import Path

from prefect import flow  # type: ignore
from gb_prefect.tasks.is_reference import is_reference  # pylint: disable=wrong-import-position


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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--file-path", required=True, help="Path to the input file to check.")
    parser.add_argument("--output-path", required=True, help="Directory where output and logs will be saved.")
    parser.add_argument("--enscode", required=True, help="Path to ENSCODE directory.")
    parser.add_argument(
        "--asm_venv", required=True, help="Path to the assembly registry virtual environment."
    )
    args = parser.parse_args()

    gb_is_reference_flow(
        file_path=args.file_path,
        output_path=args.output_path,
        enscode=args.enscode,
        asm_venv=args.asm_venv,
    )
