# TAXONOMY

TAXONOMY UPDATE
This process updates the taxonomy information for a given GCA accession using the taxonomy.py script.
Inputs:
- gca: The GCA accession for which to update taxonomy information.
- attempt_update: A flag indicating whether to attempt an update (true/false).
- metadata_json: The JSON file containing the metadata information.
Outputs:
- stdout: The standard output from the taxonomy.py script.

## Process Details

| Property | Value |
|----------|-------|
| Process | `TAXONOMY` |
| Label | `'python'` |
| Tag | `${gca}` |

## Inputs

### Nextflow interface

```nextflow
tuple val(gca), val(attempt_update), path(metadata_json)
```

## Outputs

### Nextflow interface

```nextflow
stdout
```

## Source

`/Users/ftricomi/Downloads/ensembl-genes-metadata/pipelines/assembly_metadata_update/modules/taxonomy.nf`
