# Input

Like the Assembly Metadata pipeline, this pipeline doesn't take a data
file directly — it selects which already-registered GCA accessions to
check, either from the database or from a list you provide.

## Selecting which assemblies to check

Exactly one of these must be set — the pipeline fails fast with a clear
error if neither or both are provided:

| Mode | How it's selected | Relevant parameters |
|------|--------------------|----------------------|
| Screen the database | High-priority assemblies (vertebrates and specific BioProjects) released after a given date | `--screen_date YYYY-MM-DD` |
| Custom list | A specific set of accessions | `--gca_input --gca_list <path>` |

### Custom GCA list format

When using `--gca_input`, `--gca_list` points to a plain-text file with one
GCA accession per line:

```
GCA_000001405.29
GCA_000002035.3
GCA_000003745.2
```

Accessions must start with `GCA_`.

## Configuration files

JSON files under `data/` control database writes and, optionally, Slack
notifications. Template files with empty values are provided.

### `data/db_table_conf.json`

Maps metadata fields to database tables and write methods. Used when
backfilling taxonomy records and when registering a new taxon. Shouldn't
need to change under normal use.

### `data/slack_user.json`

Only required when `--slack_report` is set. Maps genebuilder usernames to
Slack user IDs, used to identify who to notify of a metadata change:

```json
{
    "genebuilder_username": "slack_user_id"
}
```

## Required credentials and paths

| Parameter | Description |
|-----------|-------------|
| `--metadata_params_string` | JSON string with database credentials, e.g. `{"host":"host","user":"user","password":"password","port":port,"database":"database"}` |
| `--output_dir` | Path to the directory where results will be written |

## Optional: Slack reporting

| Parameter | Description |
|-----------|-------------|
| `--slack_report` | Sends a Slack DM to the genebuilder assigned to an assembly whenever its metadata is updated |
| `--slack_user` | Path to the `data/slack_user.json` mapping (defaults to that path) |
| `--slack_params` | Slack bot connection parameters as inline JSON, e.g. `{"slack_bot_token": "xoxb-..."}` |

See [Parameters](parameters.md) for the complete list.
