# verify_stable_id_ranges.py

Verify that stable ID prefixes and numeric ranges in Ensembl core databases match the metadata registry.

## What it does

For each core database, the script:

1. Retrieves the GCA accession from the `meta` table
2. Gets the min and max `stable_id` from the `gene` table
3. Parses the prefix (letters) and numeric parts separately
4. Queries the metadata registry (`gb_assembly_metadata`) for matching `species_prefix` and `stable_space` rows
5. Compares:
   - **Prefix**: the core prefix (with trailing `G` stripped) must match the registry prefix
   - **Range**: the core min/max numbers must fall within `stable_space_start` and `stable_space_end`

## Prerequisites

- Python 3.10+
- `pymysql` or `mysql-connector-python`

## Usage

```bash
# Pass core database names directly
python verify_stable_id_ranges.py core_db_1 core_db_2 core_db_3

# Or read from a file (one database name per line)
python verify_stable_id_ranges.py --db-list cores.txt

# Write results to CSV
python verify_stable_id_ranges.py --db-list cores.txt --out results.csv
```

## Output

### Terminal

A formatted table showing core DB, GCA, prefix, min/max stable IDs, and status.

### CSV (`--out`)

| Column | Description |
| --- | --- |
| `core_db` | Core database name |
| `gca_accession` | GCA accession from the meta table |
| `prefix` | Letter prefix from the core stable IDs |
| `min_stable_id` | Full min stable ID (e.g. `ENSZSPG00000000001`) |
| `max_stable_id` | Full max stable ID |
| `min_number` | Numeric part of min stable ID |
| `max_number` | Numeric part of max stable ID |
| `status` | `OK` or description of issue(s) |

## Possible statuses

- `OK` - prefix and range match the registry
- `prefix issue` - the core prefix does not match the registry
- `range issue` - the core min or max is outside the registry range
- `no prefix row found` - no matching row in `species_prefix`
- `multiple prefix rows` - more than one row returned for this GCA
- `no stable_space row found` - no matching row in `stable_space`
- `ERROR` - exception during query (e.g. connection failure)