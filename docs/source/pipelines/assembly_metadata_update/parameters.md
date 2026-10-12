# Assembly_Metadata_Update Parameters

Automatically generated from the pipeline `nextflow_schema.json`.

## Input/output options

Define where the pipeline should save output data.

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `output_dir` | string |  | yes | Path to the output directory where results will be saved. |

## Assembly options

Options controlling which assemblies are screened for updates.

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `screen_date` | string |  | no | Custom date to retrieve assemblies released after this date. Format: YYYY-MM-DD. Required unless --gca_input is set; mutually exclusive with --gca_input. |
| `gca_input` | boolean | False | no | When true, uses the GCA list provided via --gca_list as input instead of fetching assemblies from the metadata database. Mutually exclusive with --screen_date. |
| `gca_list` | string |  | no | Path to a file containing a list of GCA accessions to process. Requires --gca_input. |

## NCBI options

Options for the NCBI Datasets API.

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `ncbi_url` | string | https://api.ncbi.nlm.nih.gov/datasets/v2 | no | Base URL for the NCBI Datasets API. |

## Database options

Options for the metadata database connection.

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `metadata_params_string` | string |  | no | JSON string with metadata database connection parameters (host, user, password, port, database). |
| `db_table_conf` | string | $projectDir/data/db_table_conf.json | no | Path to the JSON file describing how each metadata database table should be updated (method, delete key, update key). |

## Slack reporting options

Options controlling Slack notifications sent to genebuilders when an assembly's metadata changes.

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `slack_report` | boolean | False | no | When true, sends a Slack direct message to the genebuilder assigned to an assembly whenever its metadata is updated. |
| `slack_user` | string | $projectDir/data/slack_user.json | no | Path to the JSON file mapping genebuilder usernames to Slack user IDs. Requires --slack_report. |
| `slack_params` | string |  | no | JSON string with Slack bot connection parameters, e.g. '{"slack_bot_token": "xoxb-..."}'. Requires --slack_report. |
