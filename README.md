# Genebuild Metadata

`ensembl-genes-metadata` is a repository that centralizes Nextflow pipelines and other helper tools for managing metadata within Genebuild.

## Nextflow Pipelines

- [Assembly Metadata Registry](./pipelines/assembly_metadata/README.md): Nextflow pipeline to update the MySQL Genebuild assembly metadata database.
- [Assembly Metadata Updates](./pipelines/assembly_metadata_update/README.md): Nextflow pipeline to check existing assemblies against NCBI and update their metadata.

## Prefect Orchestrator

- [Genebuild Metadata Prefect](./gb_prefect/README.md): Prefect flows that submit the Nextflow pipelines above (and other scripts) as SLURM jobs.

## Web App
- Genebuild Web App: *coming soon*

## Requirements

- [Nextflow](https://www.nextflow.io/) >= 24.04.03
- Singularity >= 3.7.0
- Java 11 or later
- Access to a SLURM cluster
- Access to a MySQL instance for the Genebuild assembly metadata database
- Python 3.9+ and [Prefect](https://www.prefect.io/) (for the Prefect flows)

## Installation

```bash
git clone https://github.com/Ensembl/ensembl-genes-metadata.git
cd ensembl-genes-metadata
```

See each pipeline's own README (e.g. [assembly_metadata](./pipelines/assembly_metadata/README.md), [assembly_metadata_update](./pipelines/assembly_metadata_update/README.md)) for pipeline-specific setup and configuration.

## Contributing

Contributions are welcome. Please open an issue to discuss significant changes before submitting a pull request, and make sure any new pipeline includes its own README.

## License

This project is licensed under the [Apache License 2.0](./LICENSE)