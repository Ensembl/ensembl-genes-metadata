# SPECIES_CHECKER

SPECIES CHECKER
This process checks the species information based on the provided JSON.
Inputs:
- gca: The GCA accession.
- species_tmp: The temporary species metadata file.
- last_id: The last processed ID.
Outputs:
- gca: The GCA accession.
- ${species_tmp.baseName}.json: The updated species metadata in JSON format.

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
tuple val(gca), path(species_tmp), path(last_id)
path species_checker_script
```

## Outputs

### Nextflow interface

```nextflow
tuple val(gca), path("${species_tmp.baseName}.json")
```

## Implementation Summary

- Execute Python script

## Source

`/Users/ftricomi/Downloads/ensembl-genes-metadata/pipelines/assembly_metadata/modules/species_checker.nf`
