"""Deploy the Prefect flows running Nextflow pipelines from ENSEMBL-GENES-METADATA. 
By default, they use the genebuild-codon slurm-cli work pool. The work pool must exist 
before deploying the flows.

Every deployment has concurrency_limit 1, so two runs of the same flow never write to the 
genebuild assembly metadata DB at once (extra runs wait for the slot).

The worker clones --branch from GitHub at run time, so commit and push the branch first.

Usage: run from the repository root, with PREFECT_API_URL pointing at the Prefect server

    python gb_prefect/deployments/deploy_flows.py \\
        --branch main \\
        --outdir /path/to/output \\
        --enscode /path/to/enscode \\
        --asm-venv /path/to/asm_venv \\
        [--flows register update] [--name-suffix _test] [--dry-run]

See gb_prefect/deployments/deploy_genes_nf_flows.py for flows running ensembl-genes-nf pipelines.
"""

from gb_prefect.utils.deploy_utils import FlowDeployment, main

POOL = "genebuild-codon"

FLOWS = {
    "update": FlowDeployment(
        name="gb_metadata_update",
        entrypoint="gb_prefect/flows/gb_registry_update.py:gb_registry_update_flow",
        work_pool=POOL,
        description="Assembly metadata update",
        tags=["metadata_update"],
        concurrency_limit=1,
    ),
    "by_date": FlowDeployment(
        name="gb_metadata_update_per_date",
        entrypoint="gb_prefect/flows/gb_registry_update_by_date.py:gb_registry_update_by_date_flow",
        work_pool=POOL,
        description="Assembly metadata by_date",
        tags=["metadata_update"],
        concurrency_limit=1,
    ),
    "register": FlowDeployment(
        name="gb_metadata_registry",
        entrypoint="gb_prefect/flows/gb_registry.py:gb_registry_flow",
        work_pool=POOL,
        description="Assembly metadata register",
        tags=["metadata_registry"],
        concurrency_limit=1,
    ),
    "register_by_gca": FlowDeployment(
        name="gb_metadata_registry_per_gca",
        entrypoint="gb_prefect/flows/gb_registry_by_gca.py:gb_registry_by_gca_flow",
        work_pool=POOL,
        description="Assembly metadata register_by_gca",
        tags=["metadata_registry"],
        concurrency_limit=1,
    ),
}


if __name__ == "__main__":
    main(FLOWS, description="Deploy the assembly registry/update flows (ensembl-genes-metadata pipelines).")
