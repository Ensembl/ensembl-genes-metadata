# Output

The pipeline's primary effect is updating existing records in place in the
Genebuild assembly metadata database whenever a discrepancy against NCBI is
found. Alongside that, it writes report files to `<output_dir>`.

## Final outputs

Located directly in `<output_dir>/`:

| File | Description |
|------|-------------|
| `deleted_GCAS_to_add.csv` | GCA accessions deleted from the metadata database due to missing records and no active genebuild annotation |
| `to_manually_check_GCAS.csv` | GCA accessions with missing records that could not be auto-deleted (an active genebuild annotation exists) and require manual review |
| `report_track.csv` | All metadata comparison results, one row per check performed, with header `assembly,check_type,reporting,previous_value,new_value` — covering `asm_status`, `refseq_check`, `asm_metrics`, `asm_name_check`, `bioproject_check`, `taxon_id_check`, `scientific_name_check` and `common_name_check` |

## Intermediate outputs per accession

Located in `<output_dir>/nextflow_output/<GCA>/`:

| File | Description |
|------|-------------|
| `<GCA>_metadata.json` | Assembly report retrieved from the NCBI Datasets API |
| `taxonomy_<taxon_id>.json` | Taxonomy lineage retrieved from the NCBI Taxonomy API (produced when backfilling taxonomy records or registering a new taxon) |
| `taxonomy_<taxon_id>.last_id` | Database ID assigned to the taxonomy/taxonomy_name record |
| `slack_reporting.log` | Log of the Slack notification sent for the accession (only when `--slack_report` is set) |

## Reading `report_track.csv`

Each row records one comparison outcome: the accession checked, which
field was compared (`check_type`), whether it was actually changed in the
database (`reporting`), and the previous/new values. This is the file to
consult to see exactly what changed in a given run.
