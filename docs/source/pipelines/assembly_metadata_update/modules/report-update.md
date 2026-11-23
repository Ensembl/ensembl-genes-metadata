# REPORT_UPDATE

REPORT UPDATE
This process updates the report for a given GCA accession using the report_update.py script.
Inputs:
- gca: The GCA accession for which to update the report.
- check: The type of check being performed.
- reporting: A flag indicating whether reporting is enabled (true/false).
- old_value: The previous value before the update.
- new_value: The new value after the update.
Outputs:
- slack_reporting.log: The log file containing the Slack reporting information.

## Process Details

| Property | Value |
|----------|-------|
| Process | `REPORT_UPDATE` |
| Label | `'python'` |
| Tag | `${gca}` |
| Publish directory | `"${params.output_dir}/nextflow_output/${gca}", mode: 'copy'` |

## Inputs

### Nextflow interface

```nextflow
tuple val(gca), val(check), val(reporting), val(old_value), val(new_value)
```

## Outputs

### Nextflow interface

```nextflow
path "slack_reporting.log", emit: asm_file
```

## Source

`/Users/ftricomi/Downloads/ensembl-genes-metadata/pipelines/assembly_metadata_update/modules/report_update.nf`
