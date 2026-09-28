# Input

The pipeline does not take a file of records to process by default — it
discovers which GCA accessions to register itself, either by querying NCBI
or from a list you provide.

## Selecting which assemblies to fetch

Exactly one strategy is used per run:

| Mode | How it's selected | Relevant parameters |
|------|--------------------|----------------------|
| Automatic (default) | Assemblies released since the last update date recorded in the database | none — just run the pipeline |
| Custom date | Assemblies released after a given date | `--date MM/DD/YYYY` |
| Full scan | All assemblies released since 01/01/2019 | `--full_screen` |
| Custom list | A specific set of accessions, skipping the NCBI query entirely | `--add_gca --gca_list <path>` |

`--date` and `--full_screen` are mutually exclusive. See
[Parameters](parameters.md) for the full list.

### Custom GCA list format

When using `--add_gca`, `--gca_list` points to a plain-text file with one
GCA accession per line:

```
GCA_000001405.29
GCA_000002035.3
GCA_000003745.2
```

Accessions must start with `GCA_`. Already-registered assemblies are
filtered out automatically regardless of which selection mode is used.

## Configuration files

Two JSON files under `data/` control how the pipeline talks to NCBI (`data/ncbi_params.json`) and the
database (`data/db_table_conf.json`). Template files with default values are provided; both should
rarely need changing.

### `data/ncbi_params.json`

Filters passed to the NCBI Datasets API — which assembly sources, levels
and annotation status to include:

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

Maps assembly metadata fields to the database tables and write methods
used to store them.

## Required credentials and paths

| Parameter | Description |
|-----------|-------------|
| `--metadata_params_string` | JSON string with database credentials, e.g. `{"host":"host","user":"user","password":"password","port":port,"database":"database"}` |
| `--enscode` | Path to the directory containing Ensembl repositories (`$ENSCODE`) — required so the pipeline can locate `ensembl-genes` |
| `--output_dir` | Path to the directory where results will be written |
