# gb_metadata

Shared Python package used by the Genebuild metadata Nextflow pipelines, plus a couple of standalone
command-line tools. Installed as the importable package **`gb_metadata`** — `pyproject.toml` maps the
package name `gb_metadata` to this directory (`src/python`).

## Modules

| Module | Purpose |
|---|---|
| `db_utils.py` | Shared MySQL connection helpers (`execute_query`, `execute_write`, `fetch_one_row`) used by the pipeline `bin/` scripts in both [`pipelines/assembly_metadata/`](../../pipelines/assembly_metadata/) and [`pipelines/assembly_metadata_update/`](../../pipelines/assembly_metadata_update/). |
| `utils.py` | `connection_api()` — generic HTTP GET with retry/backoff, used for NCBI and ToLID API calls. |
| `species_checker.py` | Adds `species_taxon_id`/`parlance_name` to a species record and retrieves taxonomy classification/names for the `taxonomy`/`taxonomy_name` tables from the NCBI taxonomy API. Used both as a pipeline module (`--json-path`, `--ncbi_url`, `--enscode`) and as a standalone taxonomy-backfill tool (`--taxonomy_update --taxon_id <id> --ncbi_url <url>`). |
| `write2db.py` | Builds and executes insert/update queries for a JSON payload against the metadata DB, driven by a table-config JSON file (`--config`). Shared by both pipelines' `WRITE2DB` processes. Deliberately kept independent from `db_utils.py` — see note below. |
| `is_reference.py` | Standalone CLI, **not used by any pipeline**: given a list of GCA accessions, queries the NCBI Datasets API to determine whether each is the reference assembly for its taxon. |

### `write2db.py` is intentionally separate

`write2db.py` has its own `execute_query`/`insert_query`/`update_query`/`retrieve_row_id` implementation
rather than using `db_utils.py`'s helpers — it's tailored to the insert/update/duplicate-row handling that
`db_utils.py`'s simpler helpers don't need to support. 

## Usage

Pipeline scripts import from this package as:

```python
from gb_metadata.db_utils import execute_query
from gb_metadata.utils import connection_api
```

Every module that talks to the metadata database takes its connection parameters (`host`, `user`,
`password`, `port`, `database`) as an inline JSON string, parsed via `type=json.loads` in argparse — not a
file path.

### Standalone tools

```bash
# Check whether a list of GCA accessions are NCBI reference assemblies
python src/python/is_reference.py GCA_004027535.1 GCA_000001405.29
python src/python/is_reference.py --file gcas.txt --output my_results.csv

# Backfill the taxonomy/taxonomy_name tables for a given taxon
python src/python/species_checker.py --taxonomy_update --taxon_id 9606 \
    --ncbi_url https://api.ncbi.nlm.nih.gov/datasets/v2
```

## Installation

This package is installed as part of the whole repo:

```bash
pip install .
```

The Nextflow pipelines run this code inside a container image, which must have `gb_metadata` installed
(`pip install .` against this repo). 
