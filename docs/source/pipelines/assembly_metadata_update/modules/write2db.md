# WRITE2DB

WRITE2DB
Shared process that writes a taxonomy JSON file to the DB via write2db.py.
Used both for the integrity-check taxonomy update and the post taxonomy-check
update path.
Inputs:
- gca: The GCA accession.
- taxonomy_json: The taxonomy JSON file to write to the DB.
- passthrough: Opaque call-site data to forward unchanged alongside the output.
Outputs:
- gca: The GCA accession.
- ${taxonomy_json.baseName}.last_id: The last processed ID.
- passthrough: The call-site data, forwarded unchanged.

## Process Details

| Property | Value |
|----------|-------|
| Process | `WRITE2DB` |
| Label | `'python'` |
| Tag | `${gca}` |
| Publish directory | `"${params.output_dir}/nextflow_output/${gca}", mode: 'copy'` |

## Inputs

### Nextflow interface

```nextflow
tuple val(gca), path(taxonomy_json), val(passthrough)
path write2db_script
```

## Outputs

### Nextflow interface

```nextflow
tuple val(gca), path("${taxonomy_json.baseName}.last_id"), val(passthrough)
```

## Implementation Summary

- Execute Python script

## Source

`/Users/ftricomi/Downloads/ensembl-genes-metadata/pipelines/assembly_metadata_update/modules/write2db.nf`
