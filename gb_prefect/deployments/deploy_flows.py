"""Deploy the Prefect flows running Nextflow pipelines, from ENSEMBL-GENES-METADATA (assembly
registry/update) and ENSEMBL-GENES-NF (statistics: BUSCO), to their slurm-cli work pools.
The work pools must exist before deploying the flows.

- registry/update flows: the metadata pool (default genebuild-metadata). Every deployment has
  concurrency_limit 1, so two runs of the same flow never write to the genebuild assembly
  metadata DB at once (extra runs wait for the slot).
- busco_gca: the busco pool (default genebuild-busco), whose own concurrency limit sizes the
  automatic dispatch. Give it its own output directory with --flow-outdir busco_gca=PATH.

Pool names are the POOLS defaults below; change one with --<key>-pool (e.g. --metadata-pool
my-test-pool), or send every selected flow to one pool with --work-pool.

The flow code always comes from this repository; an ensembl-genes-nf pipeline is found under
ENSCODE. To add a flow (e.g. OMArk), add its pool to POOLS and a FlowDeployment to FLOWS.

The worker clones --branch from GitHub at run time, so commit and push the branch first.

Usage: run from the repository root, with PREFECT_API_URL pointing at the Prefect server

    python gb_prefect/deployments/deploy_flows.py \\
        --branch main \\
        --outdir /path/to/output \\
        --flow-outdir busco_gca=/path/to/busco_output \\
        --enscode /path/to/enscode \\
        --asm-venv /path/to/asm_venv \\
        [--metadata-pool POOL] [--busco-pool POOL] [--work-pool POOL] \\
        [--flows register update busco_gca] [--name-suffix _test] [--dry-run]
"""

from gb_prefect.utils.deploy_utils import FlowDeployment, main

# Work pools by key; each key's name can be changed with --<key>-pool.
POOLS = {
    "metadata": "genebuild-metadata",
    "busco": "genebuild-busco",
    # "omark": "genebuild-omark",
}

FLOWS = {
    # ensembl-genes-metadata pipelines
    "update": FlowDeployment(
        name="gb_metadata_update",
        entrypoint="gb_prefect/flows/gb_registry_update.py:gb_registry_update_flow",
        work_pool="metadata",
        description="Assembly metadata update",
        tags=["metadata_update"],
        concurrency_limit=1,
    ),
    "by_date": FlowDeployment(
        name="gb_metadata_update_per_date",
        entrypoint="gb_prefect/flows/gb_registry_update_by_date.py:gb_registry_update_by_date_flow",
        work_pool="metadata",
        description="Assembly metadata by_date",
        tags=["metadata_update"],
        concurrency_limit=1,
    ),
    "register": FlowDeployment(
        name="gb_metadata_registry",
        entrypoint="gb_prefect/flows/gb_registry.py:gb_registry_flow",
        work_pool="metadata",
        description="Assembly metadata register",
        tags=["metadata_registry"],
        concurrency_limit=1,
    ),
    "register_by_gca": FlowDeployment(
        name="gb_metadata_registry_per_gca",
        entrypoint="gb_prefect/flows/gb_registry_by_gca.py:gb_registry_by_gca_flow",
        work_pool="metadata",
        description="Assembly metadata register_by_gca",
        tags=["metadata_registry"],
        concurrency_limit=1,
    ),
    # ensembl-genes-nf pipelines
    "busco_gca": FlowDeployment(
        name="busco-gca",
        entrypoint="gb_prefect/flows/gb_busco_gca.py:busco_gca_flow",
        work_pool="busco",
        description="BUSCO for a single GCA",
        tags=["busco"],
        cpu=2,
        memory=8,
        time_limit=4,
        separate_outdir=True,
    ),
    # "omark_db": FlowDeployment(name="omark-db", entrypoint=..., work_pool="omark", ...),
}


if __name__ == "__main__":
    main(
        FLOWS,
        POOLS,
        description="Deploy the Nextflow flows (ensembl-genes-metadata and ensembl-genes-nf pipelines).",
    )
