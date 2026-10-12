# BIOPROJECT

BIOPROJECT
This process runs the bioproject.py script to update bioproject information in the database.
Inputs:
- gca: The GCA accession number for the assembly.
- attempt_update: A flag indicating whether to attempt the update (true/false).
- metadata_json: Path to the JSON file containing metadata for the assembly.
Outputs:
- bioproject_update: The standard output from the bioproject.py script, a string indicating the result of the update operation.

## Process Details

| Property | Value |
|----------|-------|
| Process | `BIOPROJECT` |
| Label | `'python'` |
| Tag | `$gca` |

## Inputs

### Nextflow interface

```nextflow
tuple val(gca), val(attempt_update), path(metadata_json)
```

## Outputs

### Nextflow interface

```nextflow
stdout emit: bioproject_update
```

## Source

`/Users/ftricomi/Downloads/ensembl-genes-metadata/pipelines/assembly_metadata_update/modules/bioproject.nf`
