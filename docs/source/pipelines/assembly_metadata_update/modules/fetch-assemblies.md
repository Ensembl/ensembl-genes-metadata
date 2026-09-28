# FETCH_ASSEMBLIES

FETCH ASSEMBLIES
This process receive a list of GCAS accessions from the user or fetches a list from the assembly metadata database based on a screem date.
Inputs:
- screen_date: The date to filter assemblies in the database.
- gca_input: A flag indicating whether to use a user-provided list of GCA accessions (true/false).
- gca_list: Path to the file containing a list of GCA accessions (if gca_input is true).
Outputs:
- stdout: The standard output from the fetch_assemblies.py script or the grep command, which is a list of GCA accessions to be processed.

## Process Details

| Property | Value |
|----------|-------|
| Process | `FETCH_ASSEMBLIES` |
| Label | `'python'` |
| Tag | `date:$screen_date` |

## Inputs

### Nextflow interface

```nextflow
val screen_date
```

## Outputs

### Nextflow interface

```nextflow
stdout
```

## Source

`/Users/ftricomi/Downloads/ensembl-genes-metadata/pipelines/assembly_metadata_update/modules/fetch_assemblies.nf`
