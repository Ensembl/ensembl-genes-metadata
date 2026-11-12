# Assembly_Metadata Parameters

Automatically generated from the pipeline `nextflow_schema.json`.

## Input/output options

Define where the pipeline should find input data and save output data.

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `enscode` | string |  | yes | Path to the ENSCODE directory. |
| `output_dir` | string |  | yes | Path to the output directory where results will be saved. |

## Assembly options

Options controlling which assemblies are fetched.

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `taxon` | integer | 2759 | no | NCBI taxon ID used to filter assemblies. Default is 2759 (Eukaryota). |
| `date` | string |  | no | Custom date to retrieve assemblies released after this date. Format: MM/DD/YYYY. |
| `full_screen` | boolean | False | no | When true, retrieves all assemblies since 2019 regardless of date. |
| `add_gca` | boolean | False | no | When true, uses a user-provided GCA list as input instead of fetching from NCBI. |
| `gca_list` | string |  | no | Path to a file containing a list of GCA accessions to process. Requires --add_gca. |

## NCBI options

Options for the NCBI Datasets API.

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `ncbi_url` | string | https://api.ncbi.nlm.nih.gov/datasets/v2 | no | Base URL for the NCBI Datasets API. |
| `ncbi_params` | string |  | no | Path to the JSON file containing NCBI API query parameters. |

## Database options

Options for the metadata database connection.

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `metadata_params_string` | string |  | no | JSON string with metadata database connection parameters (host, user, password, port, database). |
| `db_table_conf` | string |  | no | Path to the JSON file containing database table configuration. |
| `registry_params` | string |  | no | Path to the JSON file containing Ensembl registry connection parameters. |
