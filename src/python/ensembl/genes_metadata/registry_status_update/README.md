# Registry Status Update Scripts

This directory contains Python utilities used to synchronise genebuild status information between the **new registry** and the **production database**. The scripts support identifying assemblies that require status updates and validating status changes across systems.

Recommended workflow:
1. Run `update_annotation_status.py` in test mode.
2. Apply `update_annotation_status.py` if happy.

Optional check for messy records (if records don't have proper versioning). This needs a manual check first! Only run it if you know what you are doing.
1. Run `no_version.py` in test mode to check for version mismatches.
2. Apply `no_version.py` if happy with updates.

The main user entry point is:

## `update_annotation_status.py`

This script is the primary tool for updating annotation statuses in the new registry. It performs the following tasks:

1. **Update assemblies that have been handed-over or released**

   * Checks assemblies present in production metadata DB.
   * Updates statuses, finds and fills missing information in the new registry if any.

The output includes a summary of updates applied and logs for missing or inconsistent fields.


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

1. Fetch current statuses from the new registry.
2. Compare with production DB.
3. Detect status changes and missing metadata.
4. Apply updates to the new registry.
5. Log results and produce a final summary.


## Arguments (`update_annotation_status.py`)

The script accepts the following command-line arguments:

| Argument                | Type                  | Description                                                                                                                                                                                                 |
|-------------------------| --------------------- |-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `-p`, `--password`      | required              | MySQL password for the write user on the genebuild registry.                                                                                                                                                |
| `-t`, `--test`          | flag                  | Runs in test mode. No updates are written to the registry; all intended updates are printed.                                                                                                                |

Notes:

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
