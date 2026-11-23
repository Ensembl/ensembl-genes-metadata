# INTEGRITY_CHECKER

INTEGRITY CHECKER
This process checks the integrity of the metadata for a given GCA accession using the integrity_checker.py script.
Inputs:
- gca: The GCA accession for which to check integrity.
Outputs:
- stdout: The standard output from the integrity_checker.py script.

## Process Details

| Property | Value |
|----------|-------|
| Process | `INTEGRITY_CHECKER` |
| Label | `'python'` |
| Tag | `${gca}` |

## Inputs

### Nextflow interface

```nextflow
val gca
```

## Outputs

### Nextflow interface

```nextflow
tuple val(gca), stdout
```

## Source

`/Users/ftricomi/Downloads/ensembl-genes-metadata/pipelines/assembly_metadata_update/modules/integrity_checker.nf`
