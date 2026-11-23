# FETCH_GCA

FETCH_GCA
This process fetches the GCA accessions for a given taxon from the NCBI
Inputs:
- taxon: The taxon ID for which to fetch GCA accessions.
- last_update: The date of the last update to filter new assemblies.
Outputs:
- stdout: The standard output containing the fetched GCA accessions.

## Process Details

| Property | Value |
|----------|-------|
| Process | `FETCH_GCA` |
| Label | `'python'` |
| Tag | `update:${last_update}` |
| Publish directory | `"${params.output_dir}/nextflow_output/", mode: 'copy'` |

## Inputs

### Nextflow interface

```nextflow
val taxon
val last_update
```

## Outputs

### Nextflow interface

```nextflow
path "assemblies_to_register.txt", emit: asm_file
stdout emit: gca
```

## Source

`/Users/ftricomi/Downloads/ensembl-genes-metadata/pipelines/assembly_metadata/modules/fetch_gca.nf`
