# FETCH_METADATA

FETCH METADATA
This process fetches metadata for a given GCA accession using the fetch_metadata.py script.
Inputs:
- gca: The GCA accession for which to fetch metadata.
Outputs:
- stdout: The standard output from the fetch_metadata.py script.

## Process Details

| Property | Value |
|----------|-------|
| Process | `FETCH_METADATA` |
| Label | `'python'` |
| Tag | `${gca}` |
| Publish directory | `"${params.output_dir}/nextflow_output/${gca}", mode: 'copy'` |

## Inputs

### Nextflow interface

```nextflow
val gca
```

## Outputs

### Nextflow interface

```nextflow
tuple val(gca), stdout, path("${gca}_metadata.json"), emit: metadata_json
```

## Source

`/Users/ftricomi/Downloads/ensembl-genes-metadata/pipelines/assembly_metadata_update/modules/fetch_metadata.nf`
