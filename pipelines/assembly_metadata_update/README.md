# Genebuild Assembly Metadata Update Pipeline

The pipeline assesses the integrity of records associated with a pre-selected list of assembly accessions. It flags corrupted records and cases where taxonomy-related information is missing or incomplete. Additionally, it connects to the Assembly Metadata Database and the NCBI API to compare metadata and detect any changes. Relevant fields are updated accordingly, and users may be notified via Slack in specific cases.

---

## Table of Contents

- [Requirements](#requirements)
- [Pipeline Overview](#pipeline-overview)
- [Configuration Files](#configuration-files)
- [Parameters](#parameters)
- [Running the Pipeline](#running-the-pipeline)
- [Outputs](#outputs)

---

## Requirements

- [Nextflow](https://www.nextflow.io/) >= 24.04.03
- Singularity >= 3.7.0
- Access to a SLURM cluster
- Access to the Ensembl assembly metadata MySQL database
- A Slack bot token and app, if using `--slack_report`

---

## Pipeline Overview

The pipeline runs in two phases for every GCA accession: an **integrity check** against the metadata database, followed by a **metadata comparison and update** against NCBI.

### Phase 1 — Integrity check

| Step | Module | Description |
|------|--------|-------------|
| 1 | `FETCH_ASSEMBLIES` | Builds the list of GCAs to check: from a user-provided list (`--gca_input`/`--gca_list`), or from the metadata database (high-priority assemblies — vertebrates and specific BioProjects — released after `--screen_date`) |
| 2 | `INTEGRITY_CHECKER` | For each accession, checks that the expected rows exist in `assembly_metrics`, `organism`, `bioproject`, `taxonomy`, `species` and `taxonomy_name`, and classifies it as `correct`, `taxonomy_update`, `delete`, or `check` (see below) |
| 3 | `INTEGRITY_TAXONOMY` | For `taxonomy_update` accessions, fetches the NCBI taxonomy lineage for the assembly's taxon ID |
| 4 | `INTEGRITY_WRITE2DB` | Writes the missing `taxonomy`/`taxonomy_name` rows to the metadata database |

`INTEGRITY_CHECKER` classifies each accession as:

| Status | Meaning | Outcome |
|--------|---------|---------|
| `correct` | All expected records exist | Proceeds to Phase 2 |
| `taxonomy_update` | All records exist except `taxonomy_name` | Backfilled via `INTEGRITY_TAXONOMY`/`INTEGRITY_WRITE2DB`, then proceeds to Phase 2 |
| `delete` | Records are missing and there is no active genebuild annotation | Deleted from the database and listed in `deleted_GCAS_to_add.csv` |
| `check` | Records are missing but an active genebuild annotation exists | Listed in `to_manually_check_GCAS.csv` for manual review; not processed further |

### Phase 2 — Metadata comparison and update

| Step | Module | Description |
|------|--------|-------------|
| 5 | `FETCH_METADATA` | Retrieves the latest assembly report for the accession from the NCBI Datasets API |
| 6 | `ASSEMBLY_STATUS` | Compares assembly status (current/suppressed/replaced) between the database and NCBI |
| 6 | `ASSEMBLY_REFSEQ` | Compares the RefSeq accession/paired assembly between the database and NCBI |
| 6 | `ASSEMBLY_METRICS` | Compares assembly metrics (e.g. contig N50, scaffold count) between the database and NCBI |
| 6 | `ASSEMBLY_NAME` | Compares the assembly name between the database and NCBI |
| 6 | `BIOPROJECT` | Compares the associated BioProject between the database and NCBI |
| 7 | `TAXONOMY_CHECK` | Compares the assembly's taxon ID against NCBI; routes to step 8 if it matches (or the new taxon already exists in `species`), otherwise to step 9 |
| 8 | `TAXONOMY` | Updates taxon ID, scientific name and common name directly |
| 9 | `SPECIES_CHECKER` | Fetches the NCBI taxonomy lineage for a new taxon ID not yet in the registry |
| 10 | `WRITE2DB` | Writes the new taxonomy/species records to the database |
| 11 | `TAXONOMY` (as `NEW_TAXONOMY`) | Updates taxon ID, scientific name and common name for the newly-registered taxon |
| 12 | `REPORT_UPDATE` | When `--slack_report` is set, sends a Slack DM to the genebuilder assigned to the assembly for `asm_status` and `refseq_check` updates |

Steps 6 run in parallel for every accession that reached Phase 2. Each comparison step updates the database in place when a difference is found and emits a line (`accession, check_type, reporting, previous_value, new_value`) used to build the run report.

![Pipeline diagram](./assembly_metadata_update_animated.svg)

---

## Configuration Files

The pipeline requires JSON configuration files located in the `data/` directory. Template files with empty values are provided.


### `data/db_table_conf.json`

Maps metadata fields to the corresponding database tables and write methods. Used when backfilling taxonomy records (`INTEGRITY_WRITE2DB`) and when registering a new taxon (`WRITE2DB`). This file should not need to be modified under normal use.

### `data/slack_user.json`

Maps genebuilder usernames to Slack user IDs, used by `REPORT_UPDATE` to identify who to notify. Only required when `--slack_report` is set.

```json
{
    "genebuilder_username": "slack_user_id"
}
```

---

## Parameters

### Required

| Parameter | Description |
|-----------|-------------|
| `--output_dir` | Path to the directory where pipeline results will be stored |
| `--metadata_params_string` | String-json with database credentials file |


### Default options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--ncbi_url` | `https://api.ncbi.nlm.nih.gov/datasets/v2` | NCBI Datasets API base URL |
| `--db_table_conf` | `data/db_table_conf.json` | Path to database table configuration file |

### Assembly options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--screen_date` | `null` | Retrieve assemblies released after this date. Format: `YYYY-MM-DD`. Mutually exclusive with `--gca_input` |
| `--gca_input` | `false` | When set, uses the GCA list provided via `--gca_list` as input instead of fetching assemblies from the metadata database. Requires `--gca_list`. Mutually exclusive with `--screen_date` |
| `--gca_list` | `null` | Path to a plain-text file of GCA accessions to check (one per line, must start with `GCA_`). Requires `--gca_input` |

> **Note:** Exactly one of `--screen_date` or `--gca_input`/`--gca_list` must be provided — the pipeline fails fast with a clear error if neither or both are set.


### Slack reporting options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--slack_report` | `false` | When set, sends a Slack DM to the genebuilder assigned to an assembly whenever its metadata is updated. Requires `--slack_user` and `--slack_params` |
| `--slack_user` | `data/slack_user.json` | Path to the JSON file mapping genebuilder usernames to Slack user IDs |
| `--slack_params` | `null` | Slack bot connection parameters, as an inline JSON string, e.g. `'{"slack_bot_token": "xoxb-..."}'` |

---

## Running the Pipeline

### Get help

```bash
nextflow run $ENSCODE/ensembl-genes-metadata/pipelines/assembly_metadata_update/main.nf --help
```

### Screen the database since a given date

```bash
nextflow run $ENSCODE/ensembl-genes-metadata/pipelines/assembly_metadata_update/main.nf \
    --metadata_params_string '{"host":"host","user":"user", "password":"password", "port": port, "database" : "database"}' \
    --output_dir /path/to/output/dir \
    --screen_date 2024-01-01
```

### Custom GCA list — check specific accessions

Provide a plain-text file with one GCA accession per line:

```
GCA_000001405.29
GCA_000002035.3
GCA_000003745.2
```

```bash
nextflow run $ENSCODE/ensembl-genes-metadata/pipelines/assembly_metadata_update/main.nf \
    --metadata_params_string '{"host":"host","user":"user", "password":"password", "port": port, "database" : "database"}' \
    --output_dir /path/to/output/dir \
    --gca_input \
    --gca_list /path/to/gca_list.txt 
```

### With Slack reporting enabled

```bash
nextflow run $ENSCODE/ensembl-genes-metadata/pipelines/assembly_metadata_update/main.nf \
    --metadata_params_string '{"host":"host","user":"user", "password":"password", "port": port, "database" : "database"}' \
    --output_dir /path/to/output/dir \
    --screen_date 2024-01-01 \
    --slack_report \
    --slack_params '{"slack_bot_token": "xoxb-your-token"}' \
    --metadata_params_string '{"host":"your-db-host","user":"your-db-user","password":"your-db-password","port":3306,"database":"your-db-name"}'
```

---

## Outputs

### Final outputs

Located directly in `<output_dir>/`:

| File | Description |
|------|-------------|
| `deleted_GCAS_to_add.csv` | GCA accessions deleted from the metadata database due to missing records and no active genebuild annotation |
| `to_manually_check_GCAS.csv` | GCA accessions with missing records that could not be auto-deleted (an active genebuild annotation exists) and require manual review |
| `report_track.csv` | All metadata comparison results, with header `assembly,check_type,reporting,previous_value,new_value` — one row per check performed (`asm_status`, `refseq_check`, `asm_metrics`, `asm_name_check`, `bioproject_check`, `taxon_id_check`, `scientific_name_check`, `common_name_check`) |

### Intermediate outputs per accession

Located in `<output_dir>/nextflow_output/<GCA>/`:

| File | Description |
|------|-------------|
| `<GCA>_metadata.json` | Assembly report retrieved from the NCBI Datasets API |
| `taxonomy_<taxon_id>.json` | Taxonomy lineage retrieved from the NCBI Taxonomy API (produced when backfilling taxonomy records or registering a new taxon) |
| `taxonomy_<taxon_id>.last_id` | Database ID assigned to the taxonomy/taxonomy_name record |
| `slack_reporting.log` | Log of the Slack notification sent for the accession (only when `--slack_report` is set) |
