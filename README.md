# Genebuild Metadata

`ensembl-genes-metadata` is a repository that centralizes Nextflow pipelines and other helper tools for managing metadata within Genebuild.

## Nextflow Pipelines

- [Assembly Metadata Registry](./assembly_metadata/README.md): Nextflow pipeline to update the MySQL Genebuild assembly metadata database.
- Assembly Metadata Updates: *coming soon*

## Prefect

- **Genebuild Metadata Prefect**: *coming soon* — Prefect flows for genebuild metadata management.


## Requirements

- [Nextflow](https://www.nextflow.io/) (>= 22.x)
- Java 11 or later
- Access to a MySQL instance for the Genebuild assembly metadata database
- Python 3.9+ and [Prefect](https://www.prefect.io/) (for the Prefect flows)

## Installation

```bash
git clone https://github.com/Ensembl/ensembl-genes-metadata.git
cd ensembl-genes-metadata
```

See each pipeline's own README (e.g. [assembly_metadata](./assembly_metadata/README.md)) for pipeline-specific setup and configuration.

## Contributing

Contributions are welcome. Please open an issue to discuss significant changes before submitting a pull request, and make sure any new pipeline includes its own README.

## License

This project is licensed under the [Apache License 2.0](./LICENSE)