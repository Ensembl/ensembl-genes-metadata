# SPECIES_CHECKER

SPECIES CHECKER
This process checks the species information for a given GCA accession using the species_checker.py script.
Inputs:
- gca: The GCA accession for which to check species information.
- attempt_update: A flag indicating whether to attempt an update (true/false).
- metadata_json: The JSON file containing the metadata information.
- old_taxon_id: The previous taxon ID before the update.
- new_taxon_id: The new taxon ID after the update.
Outputs:
- taxonomy_${new_taxon_id}.json: The JSON file containing the updated taxonomy information.

## Process Details

| Property | Value |
|----------|-------|
| Process | `SPECIES_CHECKER` |
| Label | `'python'` |
| Tag | `${gca}` |
| Publish directory | `"${params.output_dir}/nextflow_output/${gca}", mode: 'copy'` |

## Inputs

### Nextflow interface

```nextflow
tuple val(gca), val(attempt_update), path(metadata_json), val(old_taxon_id), val(new_taxon_id)
path species_checker_script
```

## Outputs

### Nextflow interface

```nextflow
tuple val(gca), val(attempt_update), path(metadata_json), path("taxonomy_${new_taxon_id}.json")
```

## Implementation Summary

- Execute Python script

## Source

`/Users/ftricomi/Downloads/ensembl-genes-metadata/pipelines/assembly_metadata_update/modules/species_checker.nf`
