# Genebuild Assembly Metadata Update Pipeline

The pipeline assesses the integrity of records associated with a pre-selected list of assembly accessions. It flags corrupted records and cases where taxonomy-related information is missing or incomplete. Additionally, it connects to the Assembly Metadata Database and the NCBI API to compare metadata and detect any changes. Relevant fields are updated accordingly, and users may be notified via Slack in specific cases.


## Requirements 
- [Nextflow](https://www.nextflow.io/) >= 24.04.03
- Singularity >= 3.7.0
- Access to a SLURM cluster
- Access to the Ensembl assembly metadata MySQL database

## Pipeline Overview




## Configuration Files




## Parameters




## Running the Pipeline

### Custome GCA list - check specific accessions

Provide a plain-text file with one GCA accession per line:

```
GCA_000001405.29
GCA_000002035.3
GCA_000003745.2
```

```bash
nextflow -C $ENSCODE/ensembl-genes-metadata/pipelines/assembly_metadata_update/nexflow.config \
    run $ENSCODE/ensembl-genes-metadata/pipelines/assembly_metadata_update/main.nf \
    --output_dir /path/to/output/dir \
    --gca_list /path/to/gca_list.txt \
    --gca_input true \
    --slack_report true
```

## Outputs







