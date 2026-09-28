"""Deploys the `slurm-cli` smoke-test flow to a work pool, for one-off validation.

Example:

    python gb_prefect/worker/deploy_smoke_test.py \
        --working-dir /nfs/production/flicek/ensembl/genebuild/<you>/slurm_cli_smoke_test \
        --asm-venv /path/to/asm_venv

Not a production deployment script -- delete once the slurm-cli worker is proven
end to end against a real Slurm job.
"""

import argparse

from gb_prefect.worker.smoke_test_flow import slurm_cli_smoke_test_flow

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pool", default="codon-slurm-smoke-pool", help="Work pool to deploy to (type slurm-cli)."
    )
    parser.add_argument(
        "--working-dir",
        required=True,
        help="Directory for the sbatch script and Slurm output/error files.",
    )
    parser.add_argument(
        "--asm-venv", required=True, help="Path to the venv with gb_prefect + prefect installed."
    )
    parser.add_argument("--cpu", type=int, default=1)
    parser.add_argument("--memory", type=int, default=1, help="Memory in GB.")
    parser.add_argument("--time-limit", type=int, default=1, help="Wall time in hours.")
    parser.add_argument("--partition", default=None, help="Slurm partition, optional.")
    args = parser.parse_args()

    slurm_cli_smoke_test_flow.deploy(
        name="slurm-cli-smoke-test",
        work_pool_name=args.pool,
        job_variables={
            "cpu": args.cpu,
            "memory": args.memory,
            "time_limit": args.time_limit,
            "partition": args.partition,
            "working_dir": args.working_dir,
            "setup_commands": [f"source {args.asm_venv}/bin/activate"],
        },
    )
