# GET_TOLID

GET_TOLID
This process retrieves the TOLID for a given GCA accession from id.tol.sanger.ac.uk
Inputs:
- gca: The GCA accession for which to fetch the TOLID.
Outputs:
- gca_tolid.json: A JSON file containing the TOLID information for the given GCA accession.

## Process Details

| Property | Value |
|----------|-------|
| Process | `GET_TOLID` |
| Label | `'python'` |
| Tag | `$gca` |
| Publish directory | `"${params.output_dir}/nextflow_output/$gca", mode: 'copy'` |

## Inputs

### Nextflow interface

```nextflow
tuple val(gca), path(last_id)
```

## Outputs

### Nextflow interface

```nextflow
tuple val(gca), path("${gca}_tolid.json")
```

## Source

`/Users/ftricomi/Downloads/ensembl-genes-metadata/pipelines/assembly_metadata/modules/get_tolid.nf`
