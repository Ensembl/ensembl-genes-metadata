# ASSEMBLY_NAME

ASSEMBLY NAME
This process runs the assembly_name.py script to update assembly names in the database.
Inputs:
- gca: The GCA accession number for the assembly.
- attempt_update: A flag indicating whether to attempt the update (true/false).
- metadata_json: Path to the JSON file containing metadata for the assembly.
Outputs:
- asm_name_update: The standard output from the assembly_name.py script, a string indicating the result of the update operation.

## Process Details

| Property | Value |
|----------|-------|
| Process | `ASSEMBLY_NAME` |
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
stdout emit: asm_name_update
```

## Source

`/Users/ftricomi/Downloads/ensembl-genes-metadata/pipelines/assembly_metadata_update/modules/assembly_name.nf`
