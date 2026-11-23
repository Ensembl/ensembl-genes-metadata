# REPORT

REPORT
This process generates a report based on the provided GCA list and last update date.
Inputs:
- gca_list: The list of GCA accessions to include in the report.
- last_update: The date of the last update to filter new assemblies.
Outputs:
- report.txt: The generated report file.
- gca_to_run_ncbi.csv: A CSV file containing formatted to be input of the BUSCO Nextflow pipeline

## Process Details

| Property | Value |
|----------|-------|
| Process | `REPORT` |
| Label | `'python'` |
| Publish directory | `"${params.output_dir}/nextflow_output/", mode: 'copy'` |

## Inputs

### Nextflow interface

```nextflow
path gca_list
val last_update
```

## Outputs

### Nextflow interface

```nextflow
path "report.txt"
path "gca_to_run_ncbi.csv"
```

## Source

`/Users/ftricomi/Downloads/ensembl-genes-metadata/pipelines/assembly_metadata/modules/report.nf`
