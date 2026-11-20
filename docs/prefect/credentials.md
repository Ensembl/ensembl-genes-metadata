# Credentials

Credentials (database passwords, Slack tokens) never appear in deployment parameters, command lines,
logs or artifacts. They live in Prefect **Secret blocks** on the server and are loaded by the flow at
run time.

## The Secret blocks

| Block | Holds | Used by |
|---|---|---|
| `gb-metadata-db-params` | the assembly metadata DB connection (JSON: `host`, `port`, `user`, `password`, `database`) | every flow that reads or writes the metadata DB |
| `gb-slack-params` | the Slack bot connection (JSON) | the update flows, when `slack_report` is on |
| `nf-metro-token` | the nf-metro server's token | live maps ([Operations](operations.md)) |

Create or replace them with
[`create_secrets.py`](https://github.com/Ensembl/ensembl-genes-metadata/blob/main/gb_prefect/setup/create_secrets.py)
(its docstring shows the JSON formats); the token is prompted for, not passed on the command line:

```bash
python gb_prefect/setup/create_secrets.py --metadata-json ~/secrets/metadata_params.json --slack-json ~/secrets/slack_params.json
```

```bash
python gb_prefect/setup/create_secrets.py --metro-token
```

Add `--overwrite` to replace existing blocks. The JSON files aren't needed afterwards.

## In a flow

```python
from gb_prefect.utils.credentials_utils import PipelineCredentials, resolve_credentials

credentials = resolve_credentials(credentials, dry_run=dry_run, metadata_secret_block=metadata_secret_block)
credentials.metadata_params_string   # the DB connection as a JSON string: for a pipeline or script
credentials.metadata_db_params       # the same as a dict: for code connecting to the DB itself
```

- Give the flow `credentials: Optional[PipelineCredentials] = None` and
  `metadata_secret_block: Optional[str] = None` parameters. Deployed runs leave both unset, so the
  standard block is loaded; someone calling the flow by hand can pass
  `credentials=PipelineCredentials(metadata_params_string='{...}')` and skip the server.
- `slack_report=True` also loads the Slack block (`credentials.slack_params`).
- `dry_run=True` loads nothing and returns placeholders, so a dry run needs no server.
- Load credentials *inside* any `try` that handles failures, so a missing block is treated like any
  other failure (see `busco_gca_flow`).

## Passing them on

| To | How | Where it ends up |
|---|---|---|
| a Nextflow pipeline | `NextflowRun.secret_params={"metadata_params_string": ...}` | an owner-only (`0600`) `-params-file`; redacted from the command record, logs and artifact |
| a script | `ScriptRun.secret_env={"METADATA_PARAMS": ...}`, read by the script from its environment | the script's environment only; redacted from logs and artifact, not written to the command record |
| your own Python | `credentials.metadata_db_params` | in memory |

Never pass credentials in `params` or `argv`: command lines are visible in process lists and are
written to Nextflow's `.nextflow.log`.

!!! warning "Inside the pipeline"
    The params file keeps credentials out of what `gb_prefect` writes, but a pipeline that puts a
    parameter into a process's script (e.g. `'${params.metadata_params_string}'`) still writes it to
    that process's `work/*/.command.sh`. Avoiding that is a change in the pipeline.

## Adding a new kind of credential

For example, read-write credentials for a core database:

1. Choose a block name, e.g. `gb-core-db-params`, and add a constant for it next to the others in
   [`utils/credentials_utils.py`](https://github.com/Ensembl/ensembl-genes-metadata/blob/main/gb_prefect/utils/credentials_utils.py).
2. Add an option to `create_secrets.py` to store it (from a JSON file, or prompted with `getpass`
   for a single value; never as a command-line value).
3. In the flow, load it with `Secret.load(<block>).get()` (from `prefect.blocks.system`), or extend
   `resolve_credentials` if several flows need it; skip loading on a dry run.
4. Pass it on as above (`secret_params` or `secret_env`).
5. List the block in the table above and in the README's setup steps.
