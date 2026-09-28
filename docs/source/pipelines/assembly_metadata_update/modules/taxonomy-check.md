# TAXONOMY_CHECK

TAXONOMY CHECK
This process checks the taxonomy information for a given GCA accession using the taxonomy.py script.
Inputs:
- gca: The GCA accession for which to check taxonomy information.
- attempt_update: A flag indicating whether to attempt an update (true/false).
- metadata_json: The JSON file containing the metadata information.
Outputs:
- OLD_TAXON_ID: The previous taxon ID before the update.
- NEW_TAXON_ID: The new taxon ID after the update.
- STATUS: The status of the taxonomy check.

## Process Details

| Property | Value |
|----------|-------|
| Process | `TAXONOMY_CHECK` |
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
tuple val(gca), val(attempt_update), path(metadata_json), env('OLD_TAXON_ID'), env('NEW_TAXON_ID'), env('STATUS'), emit: taxonomy_check
```

## Source

`/Users/ftricomi/Downloads/ensembl-genes-metadata/pipelines/assembly_metadata_update/modules/taxonomy_check.nf`
