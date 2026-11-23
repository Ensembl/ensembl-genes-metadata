# WRITE2DB

WRITE2DB
Shared process that writes a single JSON/tmp file to the DB via write2db.py.
Used at each stage of the registration pipeline (assembly, metadata, species,
tolid); files that need to ride alongside for the next process in the chain
are threaded through as an opaque passthrough path list rather than being
known to this process.
Inputs:
- gca: The GCA accession.
- file_to_write: The file to write to the DB.
- passthrough: Other files (possibly empty) to forward unchanged alongside the output.
- update_flag: Whether to pass --update to write2db.py.
Outputs:
- gca: The GCA accession.
- ${file_to_write.baseName}.last_id: The last processed ID.
- passthrough: The other files, forwarded unchanged.

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
tuple val(gca), path(file_to_write), val(passthrough)
path write2db_script
val update_flag
```

## Outputs

### Nextflow interface

```nextflow
tuple val(gca), path("${file_to_write.baseName}.last_id"), val(passthrough)
```

## Implementation Summary

- Execute Python script

## Source

`/Users/ftricomi/Downloads/ensembl-genes-metadata/pipelines/assembly_metadata/modules/write2db.nf`
