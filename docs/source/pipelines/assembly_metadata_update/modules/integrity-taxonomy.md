# INTEGRITY_TAXONOMY

INTEGRITY TAXONOMY
This process checks the integrity of the taxonomy for a given GCA accession using the species_checker.py script.
Inputs:
- gca: The GCA accession for which to check taxonomy.
- taxon_id: The taxon ID associated with the GCA accession.
Outputs:
- taxonomy_${taxon_id}.json: The JSON file containing the taxonomy information.

## Process Details

| Property | Value |
|----------|-------|
| Process | `INTEGRITY_TAXONOMY` |
| Label | `'python'` |
| Tag | `${gca}` |
| Publish directory | `"${params.output_dir}/nextflow_output/${gca}", mode: 'copy'` |

## Inputs

### Nextflow interface

```nextflow
tuple val(gca), val(taxon_id)
path species_checker_script
```

## Outputs

### Nextflow interface

```nextflow
tuple val(gca), path("taxonomy_${taxon_id}.json")
```

## Implementation Summary

- Execute Python script

## Source

`/Users/ftricomi/Downloads/ensembl-genes-metadata/pipelines/assembly_metadata_update/modules/integrity_taxonomy.nf`
