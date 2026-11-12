# PARSE_METADATA

PARSE_METADATA
This process retrieves the assembly metadata for a given GCA accession from the NCBI
Inputs:
- gca: The GCA accession for which to retrieve assembly metadata.
Outputs:
- gca: The GCA accession.
- ${gca}_assembly.json: The assembly metadata in JSON format.
- ${gca}_metadata.tmp: The assembly metadata in temporary format.
- ${gca}_species.tmp: The species metadata in temporary format.

## Process Details

| Property | Value |
|----------|-------|
| Process | `PARSE_METADATA` |
| Label | `'python'` |
| Tag | `$gca` |
| Publish directory | `"${params.output_dir}/nextflow_output/$gca", mode: 'copy'` |

## Inputs

### Nextflow interface

```nextflow
val gca
```

## Outputs

### Nextflow interface

```nextflow
tuple val(gca), path("${gca}_assembly.json"), path("${gca}_metadata.tmp"), path("${gca}_species.tmp")
```

## Source

`/Users/ftricomi/Downloads/ensembl-genes-metadata/pipelines/assembly_metadata/modules/parse_metadata.nf`
