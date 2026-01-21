# Registry Status Update Scripts

This directory contains a collection of Python utilities used to synchronise genebuild status information between the **old annotation registry**, the **new registry**, and the **production database**. The scripts support identifying assemblies that require status updates, transferring missing records, and validating status changes across systems.

Recommended workflow (temporary):
1. Run `no_version.py` in test mode to check for version miss matches.
2. Apply `no_version.py` if happy with updates.
3. Run `update_annotation_status.py` in test mode.
4. Apply `update_annotation_status.py` if happy.

The main user entry point is:

## `update_annotation_status.py`

This script is the primary tool for updating annotation statuses in the new registry. It performs two core tasks:

1. **Identify assemblies missing from the new registry**

   * Retrieves assemblies from the old registry.
   * Filters those **not present** in the new registry or in the production database.
   * Copies these missing assemblies into the new registry.

2. **Update assemblies tracked in both registries**

   * Checks assemblies present in both the old and new registry.
   * Detects cases where an assembly has transitioned from an *in-progress* state to a **terminal status** and prints the accessions in the logfile

3. **Update assemblies that have been handed-over or released**

   * Checks assemblies present in production metadata DB.
   * Updates statuses, finds and fills missing information in the new registry if any.

The output includes a summary of updates applied and logs for missing or inconsistent fields.

---

## `copy_status_old_registry.py`

Retrieves the full genebuild status table from the old registry and prepares it for comparison against the new registry.

---

## `production_check.py`

Queries the production database to determine whether an assembly is already tracked there.
Provides functions for filtering out completed or already-handled assemblies.

---

## `helper.py`

Contains common helper utilities used across the scripts, including MySQL query wrappers and functions for merging and validating results.

---

## `logger_settings.py`

Defines a consistent logging setup used by all scripts.
Ensures uniform formatting of info, warning, and error messages during registry synchronisation.

---

## Workflow Summary

1. Fetch old registry statuses.
2. Compare with new registry and production DB.
3. Identify new assemblies to import.
4. Detect status changes for assemblies tracked by both registries.
5. Apply updates to the new registry.
6. Log results and produce a final summary.


## Arguments (`update_annotation_status.py`)

The script accepts the following command-line arguments:

| Argument                | Type                  | Description                                                                                                                                                                                                 |
|-------------------------| --------------------- |-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `-p`, `--password`      | required              | MySQL password for the write user on the genebuild registry.                                                                                                                                                |
| `-t`, `--test`          | flag                  | Runs in test mode. No updates are written to the registry; all intended updates are printed.                                                                                                                |
| `-am`, `--apply_method` | flag (default: False)       | Controls copying method for live and handed over assemblies from the prodiction DB to the registru. If supplied, copying **is applied**; without it, the script performs checks but does not apply changes. |
| `-or`, `--old_registry` | flag                  | When supplied, entries from the old registry are checked and optionally copied into the new registry.                                                                                                       |
| `-ao`, `--apply_old`    | flag (default: False) | Controls copying from the old registry. If supplied, copying **is applied**; without it, the script performs checks but does not apply changes.                                                             |

Notes:

* `--apply_old` and `--apply_method` behave as tests because they; the default is to **not** apply copying to the registry.
* No hostnames, ports, or database names are user inputs; these are hard-coded in the script.

---

## Example Run

### Test mode (no updates applied)

```bash
python update_annotation_status.py \
  --password mypassword \
  --test
```

### Apply production → registry updates (normal use)

```bash
python update_annotation_status.py \
  --password mypassword
```

### Include old-registry check but **do not** apply old-registry updates

```bash
python update_annotation_status.py \
  --password mypassword \
  --old_registry
```

### Include old-registry check **and apply** updates from old registry

```bash
python update_annotation_status.py \
  --password mypassword \
  --old_registry \
  --apply_old
```

### Apply method copy production → registry, but don't apply status changes (recommended usage of apply_method)

```bash
python update_annotation_status.py \
  --password mypassword \
  --apply_method \
  --test
```

