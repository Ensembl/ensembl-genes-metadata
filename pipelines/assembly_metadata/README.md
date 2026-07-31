# Genebuild Assembly Metadata Pipeline

The pipeline connects to the NCBI API to retrieve the latest eukaryotic genome assemblies and registers the associated metadata into the Ensembl assembly metadata MySQL database. Data is primarily sourced from the NCBI Datasets API, complemented by species information from the NCBI Taxonomy API and ToL IDs from the [Darwin Tree of Life API](https://id.tol.sanger.ac.uk).


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

### Ensembl dependencies

It is recommended that all repositories are cloned into the same folder (`$ENSCODE`).

| Repository | Branch | URL |
|------------|--------|-----|
| ensembl-genes | main | https://github.com/Ensembl/ensembl-genes.git |

---

## Pipeline Overview

The pipeline runs the following steps in sequence for each GCA accession:

| Step | Module | Description |
|------|--------|-------------|
| 1 | `SET_DATE` | Determines the date threshold for fetching new assemblies (from DB, custom, or full scan since 2019) |
| 2 | `FETCH_GCA` | Fetches GCA accessions from the NCBI API (or from a user-provided list) and filters out already-registered assemblies |
| 3 | `PARSE_METADATA` | Retrieves and parses assembly metadata from the NCBI Datasets API |
| 4 | `WRITE2DB_ASSEMBLY` | Writes assembly data to the `assembly` table in the metadata database |
| 5 | `UPDATE_KEYS_METADATA` | Updates foreign keys in the metadata JSON using the newly assigned DB IDs |
| 6 | `WRITE2DB_METADATA` | Writes extended metadata (metrics, bioprojects, taxonomy) to the metadata database |
| 7 | `SPECIES_CHECKER` | Validates and enriches species information using the NCBI Taxonomy API |
| 8 | `WRITE2DB_SPECIES` | Writes species data to the `species` and `organism` tables |
| 9 | `GET_TOLID` | Queries the Darwin Tree of Life API to retrieve the ToL ID for the assembly |
| 10 | `WRITE2DB_TOLID` | Writes the ToL ID to the metadata database |
| 11 | `REPORT` | Generates a summary report and a CSV file for downstream BUSCO analyses |

![Pipeline diagram](./assembly_metadata_animated.svg)

---

## Configuration Files

The pipeline requires two JSON configuration files located in the `data/` directory. Template files with empty values are provided.


### `data/ncbi_params.json`

Parameters for the NCBI Datasets API. These control which assemblies are retrieved. The defaults are suitable for most runs and should not need to be changed.

```json
{
    "filters.reference_only": "false",
    "filters.assembly_source": "genbank",
    "filters.has_annotation": "false",
    "filters.exclude_paired_reports": "false",
    "filters.exclude_atypical": "true",
    "filters.assembly_version": "current",
    "filters.assembly_level": ["scaffold", "chromosome", "complete_genome", "contig"],
    "filters.is_metagenome_derived": "metagenome_derived_exclude",
    "returned_content": "ASSM_ACC",
    "page_size": "100"
}
```

### `data/db_table_conf.json`

Maps assembly metadata fields to the corresponding database tables and write methods. This file should not need to be modified under normal use.

---

## Parameters

### Required

| Parameter | Description |
|-----------|-------------|
| `--output_dir` | Path to the directory where pipeline results will be stored |
| `--enscode` | Path to the directory containing Ensembl repositories (`$ENSCODE`) |
| `--metadata_params_string` | JSON string with database credentials |

### Assembly options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--taxon` | `2759` | NCBI taxon ID to filter assemblies. Default is 2759 (Eukaryota) |
| `--date` | `null` | Retrieve assemblies released after this date. Format: `MM/DD/YYYY`. Cannot be used together with `--full_screen` |
| `--full_screen` | `false` | When set, retrieves all assemblies since 01/01/2019. Cannot be used together with `--date` |
| `--add_gca` | `false` | When set, uses a user-provided GCA list as input instead of fetching from NCBI. Requires `--gca_list` |
| `--gca_list` | `null` | Path to a plain-text file of GCA accessions to register (one per line, must start with `GCA_`). Requires `--add_gca` |

> **Note:** `--date` and `--full_screen` are mutually exclusive. If neither is set, the pipeline uses the last regular update date stored in the metadata database (minus one day).

### Internal defaults (do not change)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--ncbi_url` | `https://api.ncbi.nlm.nih.gov/datasets/v2` | NCBI Datasets API base URL |
| `--ncbi_params` | `data/ncbi_params.json` | Path to NCBI API parameters file |
| `--db_table_conf` | `data/db_table_conf.json` | Path to database table configuration file |

---

## Running the Pipeline

### Get help

```bash
nextflow run ${ENSCODE}/ensembl-genes-metadata/pipelines/assembly_metadata/main.nf --help
```

### Default run — automatic date from database

The pipeline reads the last update date from the metadata database and fetches assemblies released since then.

```bash
nextflow -C ${ENSCODE}/ensembl-genes-metadata/pipelines/assembly_metadata/nextflow.config \
    run ${ENSCODE}/ensembl-genes-metadata/pipelines/assembly_metadata/main.nf \
    --output_dir /path/to/output/dir \
    --enscode ${ENSCODE}
```

### Manual date — fetch assemblies since a specific date

```bash
nextflow -C ${ENSCODE}/ensembl-genes-metadata/pipelines/assembly_metadata/nextflow.config \
    run ${ENSCODE}/ensembl-genes-metadata/pipelines/assembly_metadata/main.nf \
    --metadata_params_string '{"host":"host","user":"user", "password":"password", "port": port, "database" : "database"}' \
    --output_dir /path/to/output/dir \
    --enscode ${ENSCODE} \
    --date 01/01/2024
```

### Full scan — fetch all assemblies since 2019

```bash
nextflow -C ${ENSCODE}/ensembl-genes-metadata/pipelines/assembly_metadata/nextflow.config \
    run ${ENSCODE}/ensembl-genes-metadata/pipelines/assembly_metadata/main.nf \
    --metadata_params_string '{"host":"host","user":"user", "password":"password", "port": port, "database" : "database"}' \
    --output_dir /path/to/output/dir \
    --enscode ${ENSCODE} \
    --full_screen
```

### Custom GCA list — register specific accessions

Provide a plain-text file with one GCA accession per line:

```
GCA_000001405.29
GCA_000002035.3
GCA_000003745.2
```

```bash
nextflow -C ${ENSCODE}/ensembl-genes-metadata/pipelines/assembly_metadata/nextflow.config \
    run ${ENSCODE}/ensembl-genes-metadata/pipelines/assembly_metadata/main.nf \
    --metadata_params_string '{"host":"host","user":"user", "password":"password", "port": port, "database" : "database"}' \
    --output_dir /path/to/output/dir \
    --enscode ${ENSCODE} \
    --add_gca \
    --gca_list /path/to/gca_list.txt
```

---

## Outputs

### Final outputs

Located in `<output_dir>/nextflow_output/`:

| File | Description |
|------|-------------|
| `report.txt` | Summary report of newly registered assemblies: counts by assembly type/level, assemblies from relevant BioProjects, RefSeq availability, and flagged issues (invalid taxon IDs, missing BioSample or submitter info) |
| `gca_to_run_ncbi.csv` | CSV file formatted as input for the BUSCO Nextflow pipeline, listing shortlisted assemblies for genome quality assessment |
| `gca_list_to_report.txt` | List of all GCA accessions processed in the run |

### Intermediate outputs per accession

Located in `<output_dir>/nextflow_output/<GCA>/`:

| File | Description |
|------|-------------|
| `<GCA>_assembly.json` | Raw assembly metadata retrieved from NCBI |
| `<GCA>_metadata.json` | Processed metadata with updated database keys |
| `<GCA>_species.json` | Species information validated against NCBI Taxonomy |
| `<GCA>_tolid.json` | ToL ID information from the Darwin Tree of Life API |
| `<GCA>_assembly.last_id` | Database ID assigned to the assembly record |
| `<GCA>_metadata.last_id` | Database ID assigned to the metadata record |
| `<GCA>_species.last_id` | Database ID assigned to the species record |
| `<GCA>_tolid.last_id` | Database ID assigned to the ToL ID record |

