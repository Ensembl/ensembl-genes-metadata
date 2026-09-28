# UPDATE_KEYS_METADATA

UPDATE_KEYS_METADATA
This process updates the keys in the metadata JSON based on the provided temporary metadata and last processed ID
Inputs:
- gca: The GCA accession.
- metadata_tmp: The temporary metadata file.
- last_id: The last processed ID.
- species_tmp: The temporary species metadata file.
Outputs:
- gca: The GCA accession.
- ${metadata_tmp.baseName}.json: The updated metadata in JSON format.
- species_tmp: The temporary species metadata file (unchanged).

## Process Details

| Property | Value |
|----------|-------|
| Process | `UPDATE_KEYS_METADATA` |
| Label | `'python'` |
| Tag | `$gca` |
| Publish directory | `"${params.output_dir}/nextflow_output/$gca", mode: 'copy'` |

## Inputs

### Nextflow interface

```nextflow
tuple val(gca), path(metadata_tmp), path(last_id), path(species_tmp)
```

## Outputs

### Nextflow interface

```nextflow
tuple val(gca), path("${metadata_tmp.baseName}.json"), path(species_tmp)
```

## Source

`/Users/ftricomi/Downloads/ensembl-genes-metadata/pipelines/assembly_metadata/modules/update_keys_metadata.nf`
