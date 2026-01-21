"""
update_annotation_status.py

This script synchronizes the genebuild status of genome assemblies between the
Ensembl production database and the genebuild registry. It compares the registry
entries against the current production status, identifies discrepancies in status
or release date, and updates the registry accordingly.

Key Features:
- Optionally copies entries from the old registry.
- Fetches current genebuild statuses from the registry.
- Compares against production database statuses.
- Updates registry entries where the status or release date has changed.
- Supports a test mode to preview changes without applying them.

Modules:
- get_genebuild_status: Retrieves genebuild status entries from the registry.
- get_status_updates: Determines which entries require updates.
- update_genebuild_status: Applies updates to the registry.
- main: Orchestrates the workflow with optional old registry checks and test mode.

Usage:
    python update_genebuild_status.py -p <MYSQL_PASSWORD> [options]

Arguments:
    -p, --password        MySQL password for the write user.
    -t, --test            Run in test mode (no database updates applied).
    -or, --old_registry   Boolean flag to check and copy entries from the old registry.
    -sa, --stop_apply     Boolean flag; if True, do not apply updates from old registry (default: True).

Example:
    # Run in test mode
    python update_genebuild_status.py -p mypassword -t

    # Apply updates from production to registry
    python update_genebuild_status.py -p mypassword

Dependencies:
    - pandas
    - pymysql
    - copy_status_old_registry (module)
    - helper (module with mysql_fetch_data function)
    - production_check (module with check_status_production_db function)
    - logger_settings (module with get_logger function)

Notes:
- Requires network access to the Ensembl MySQL production servers.
- Outputs a 'faulty_status.csv' file for assemblies with faulty production status.
- Release dates are normalized and updated only when discrepancies are detected.
"""

import argparse
import pandas as pd
import pymysql
from copy_status_old_registry import insert_entries_from_old_registry
from helper import mysql_fetch_data
from production_check import check_status_production_db
from missing_method import add_missing_methods
from logger_settings import get_logger

logger = get_logger(__name__)

def get_genebuild_status():
    """
    Fetch genebuild status from the registry excluding 'archive' and 'live'.


    Returns:
        dict or None: First row of results as a dictionary, or None if no results.
    """
    try:
        registry_query = """
            SELECT 
                gca_accession, 
                genebuild_status_id,
                gb_status,
                genebuild_version,
                annotation_method,
                last_genebuild_update,
                date_status_update,
                release_type,
                release_date
            FROM genebuild_status
        """

        gb_status = mysql_fetch_data(
            registry_query,
            host="mysql-ens-genebuild-prod-1",
            user="ensro",
            port=4527,
            database="gb_assembly_metadata",
            password=""
        )

        logger.info(f"Found {len(gb_status)} entries in genebuild_status table.")
        gb_status = pd.DataFrame(gb_status)
        return gb_status

    except pymysql.Error as err:
        logger.error("MySQL error: %s", err)
        return


def get_status_updates(merged_df):
    """
    Determine which GCAs need their genebuild status or release date updated.
    """
    df = merged_df.copy()
    logger.info(f"Starting status update check on {len(df)} rows.")

    # Normalisation
    logger.info("Normalising status strings and dates.")
    df['gb_status'] = df['gb_status'].astype(str).str.strip().str.lower()
    df['status'] = df['status'].astype(str).str.strip().str.capitalize()

    df['release_date_registry'] = pd.to_datetime(df.get('release_date_registry', pd.NaT), errors='coerce')
    df['release_date_production'] = pd.to_datetime(df.get('release_date_production', pd.NaT), errors='coerce')
    df['last_genebuild_update_registry'] = pd.to_datetime(df.get('last_genebuild_update_registry', pd.NaT), errors='coerce')
    df['last_genebuild_update_production'] = pd.to_datetime(df.get('last_genebuild_update_production', pd.NaT), errors='coerce')
    df['date_status_update'] = pd.to_datetime(df.get('date_status_update', pd.NaT), errors='coerce')

    df['gb_status_new'] = df['gb_status']
    df['release_date_new'] = df['release_date_registry']
    df['last_genebuild_update_new'] = df['last_genebuild_update_registry']
    df['genebuild_version_new'] = df['genebuild_version']
    df['date_status_update_new'] = df['date_status_update']
    df['release_type_new'] = df['release_type']

    # 1. Faulty entries → skip
    condition_faulty = df['status'] == 'Faulty'
    faulty_count = condition_faulty.sum()
    logger.info(f"Production 'Faulty' entries: {faulty_count}")

    if faulty_count:
        faulty_path = "faulty_status.csv"
        df[condition_faulty].to_csv(faulty_path, index=False)
        logger.info(f"Saved faulty entries to: {faulty_path}")
        df = df[~condition_faulty].copy()
        logger.info(f"{len(df)} rows remain after filtering out faulty entries.")

    # 2. Production Released → make registry Live
    condition_released = (df['status'] == 'Released') & (df['gb_status'] != 'live')
    logger.info(f"Handed over→Live updates: {condition_released.sum()}")
    df.loc[condition_released, 'gb_status_new'] = 'live'
    df.loc[condition_released, 'release_date_new'] = df.loc[condition_released, 'release_date_production']
    df.loc[condition_released, 'date_status_update_new'] = pd.Timestamp.today().normalize()
    df.loc[condition_released, 'release_site_new'] = 'beta'

    # 3. Live but missing release_date
    condition_missing_date = (
        (df['gb_status'] == 'live')
        & df['release_date_registry'].isna()
        & df['release_date_production'].notna()
    )
    logger.info(f"Missing release_date to fill: {condition_missing_date.sum()}")
    df.loc[condition_missing_date, 'release_date_new'] = df.loc[condition_missing_date, 'release_date_production']
    df.loc[condition_missing_date, 'release_site_new'] = 'beta'

    # 4. Live but release_date mismatch
    condition_mismatch_date = (
        (df['gb_status'] == 'live')
        & df['release_date_registry'].notna()
        & df['release_date_production'].notna()
        & (df['release_date_registry'] != df['release_date_production'])
    )
    logger.info(f"Release_date mismatches: {condition_mismatch_date.sum()}")
    df.loc[condition_mismatch_date, 'release_date_new'] = df.loc[condition_mismatch_date, 'release_date_production']

    # 5–6. last_genebuild_update corrections
    condition_missing_update = (
        (df['gb_status'] == 'live')
        & df['last_genebuild_update_registry'].isna()
        & df['last_genebuild_update_production'].notna()
    )
    logger.info(f"Missing last_genebuild_update to fill: {condition_missing_update.sum()}")
    df.loc[condition_missing_update, 'last_genebuild_update_new'] = df.loc[condition_missing_update, 'last_genebuild_update_production']

    condition_mismatch_update = (
        (df['gb_status'] == 'live')
        & df['last_genebuild_update_registry'].notna()
        & df['last_genebuild_update_production'].notna()
        & (df['last_genebuild_update_registry'] != df['last_genebuild_update_production'])
    )
    logger.info(f"last_genebuild_update mismatches: {condition_mismatch_update.sum()}")
    df.loc[condition_mismatch_update, 'last_genebuild_update_new'] = df.loc[condition_mismatch_update, 'last_genebuild_update_production']

    # 7. Live but missing genebuild_version
    condition_missing_genebuild_version = (
            (df['gb_status'] == 'live')
            & df['genebuild_version'].isna()
            & df['gb_v_production'].notna())
    logger.info(
        f"Missing genebuild_version to fill: {condition_missing_genebuild_version.sum()}")
    df.loc[condition_missing_genebuild_version, 'genebuild_version_new'] =  df.loc[condition_missing_genebuild_version, 'gb_v_production']


    # 8. Processing/Submitted → handed_over
    condition_handed_over = df['status'].isin(['Processed', 'Processing', 'Submitted']) & (df['gb_status'] != 'handed_over')
    logger.info(f"Handed_over updates: {condition_handed_over.sum()}")
    df.loc[condition_handed_over, 'gb_status_new'] = 'handed_over'
    df.loc[condition_handed_over, 'date_status_update_new'] = pd.Timestamp.today().normalize()

    # ---- Determine changes ----
    status_changed = df['gb_status'] != df['gb_status_new']
    release_changed = ~((df['release_date_new'].isna() & df['release_date_registry'].isna()) |
                        (df['release_date_new'] == df['release_date_registry']))
    update_changed = ~((df['last_genebuild_update_new'].isna() & df['last_genebuild_update_registry'].isna()) |
                       (df['last_genebuild_update_new'] == df['last_genebuild_update_registry']))
    update_version = ~((df['genebuild_version_new'].isna() & df['genebuild_version'].isna()) |
                       (df['genebuild_version_new'] == df['genebuild_version']))



    updated_df = df[status_changed | release_changed | update_changed | update_version].copy()
    logger.info(f"Total rows requiring updates: {len(updated_df)}")

    updated_df = updated_df[
        [
            'genebuild_status_id', 'gb_status_new', 'release_date_new',
            'last_genebuild_update_new',
            'annotation_method', 'annotation_source', 'date_status_update_new', 'genebuild_version_new', 'release_site_new'
        ]
    ]

    return updated_df


def update_genebuild_status(updated_df, password):
    """
    Update the genebuild_status table with new statuses, dates info.
    Logs a summary of which columns were updated.

    Args:
        updated_df (pd.DataFrame): DataFrame with columns:
            - genebuild_status_id
            - gb_status_new
            - release_date_new (optional)
            - date_status_update_new (optional)
            - last_genebuild_update_new (optional)
        password: MySQL connection info
    """
    connection = pymysql.connect(
        host="mysql-ens-genebuild-prod-1",
        user="ensadmin",
        password=password,
        port=4527,
        database="gb_assembly_metadata",
        autocommit=True
    )

    # Track counts per column for summary
    update_counts = {
        'gb_status': 0,
        'release_date': 0,
        'date_status_update': 0,
        'last_genebuild_update': 0,
        'genebuild_version': 0,
        'release_type': 0,
    }

    try:
        with connection.cursor() as cursor:
            for _, row in updated_df.iterrows():
                genebuild_status_id = row['genebuild_status_id']
                set_clauses = []
                values = []

                # gb_status is always updated
                set_clauses.append("gb_status = %s")
                values.append(row['gb_status_new'])
                update_counts['gb_status'] += 1

                # Optional fields
                optional_fields = [
                    ('release_date_new', 'release_date'),
                    ('date_status_update_new', 'date_status_update'),
                    ('last_genebuild_update_new', 'last_genebuild_update'),
                    ('release_site_new', 'release_type'),
                    ('genebuild_version_new', 'genebuild_version'),
                ]

                updated_cols = ['gb_status']
                for df_col, db_col in optional_fields:
                    val = row.get(df_col)
                    if pd.notnull(val):
                        set_clauses.append(f"{db_col} = %s")
                        values.append(val)
                        updated_cols.append(db_col)
                        update_counts[db_col] += 1

                if not set_clauses:
                    logger.info(f"No fields to update for genebuild_status_id {genebuild_status_id}")
                    continue

                sql = f"""
                    UPDATE genebuild_status
                    SET {', '.join(set_clauses)}
                    WHERE genebuild_status_id = %s
                """
                values.append(genebuild_status_id)
                cursor.execute(sql, values)

                logger.debug(
                    f"Updated genebuild_status_id {genebuild_status_id} with fields: {', '.join(updated_cols)}"
                )

        # Summary log
        logger.info(f"Finished updating {len(updated_df)} genebuild_status rows.")
        logger.info("Update summary per column:")
        for col, count in update_counts.items():
            logger.info(f"  {col}: {count} rows updated")

    except pymysql.Error as e:
        logger.error("MySQL error: %s", e)
        raise

    finally:
        connection.close()




def main(password, test, apply_method, old_registry, apply_old):

    if old_registry:
        logger.info("Checking status in old registry.")
        gb_status = get_genebuild_status()
        copy_from_old_registry = insert_entries_from_old_registry(password, gb_status, apply_old)


    logger.info("Only using the new registry.")
    logger.info("Fetching genebuild status from new registry.")
    gb_status = get_genebuild_status()

    # Make sure we have a 'gca_accession' column
    if 'gca_accession' not in gb_status.columns:
        logger.error("No 'gca_accession' column found in registry data.")
        return

    # Convert GCAs to tuple for SQL IN clause
    gca_tuple = tuple(gb_status['gca_accession'].unique())
    if len(gca_tuple) == 1:
        gca_tuple = (gca_tuple[0],)


    logger.info(f"Getting GCA status from production DB.")
    logger.info(f"Looking for {len(gca_tuple)} accessions in production DB")
    production_status = check_status_production_db(gca_tuple)

    # Ensure consistent column names
    if 'gca_accession' not in production_status.columns:
        logger.error("Expected columns missing from production DB query.")
        return

    #Add missing method for live and handed over records
    logger.info(f"Looking for missing methods")
    missing_method = add_missing_methods(gb_status, production_status, apply_method, password)

    if production_status.empty:
        logger.warning("No production status found. Skipping merge and updates.")
        merged_df = gb_status.copy()
    else:
        logger.info("Merging dataframes.")
        production_status["gb_v_production"] = production_status["genebuild_version"]
        merge_keys = [
            "gca_accession",
            "annotation_method",
	        "genebuild_version",

        ]
        dups = gb_status.groupby(merge_keys).size()
        if (dups > 1).any():
            logger.error("Duplicate registry records detected: %s", dups[dups > 1])
            raise ValueError("Registry contains non-unique keys for GCA/method/version.")

        prod_dups = production_status.groupby(merge_keys).size()
        if (prod_dups > 1).any():
            logger.error("Duplicate production records detected: %s", prod_dups[prod_dups > 1])
            raise ValueError("Production DB contains ambiguous records.")

        logger.info("Merging registry and production dataframes.")

        merged_df = pd.merge(
            gb_status,
            production_status,
            on=merge_keys,
            how='left',
            suffixes=('_registry', '_production')
        )


    updates = get_status_updates(merged_df)

    if not test:
        update_genebuild_status(updates, password)
    else:
        logger.info("Test mode: no updates applied.")
        logger.info("Annotations that would be updated:\n%s", updates.to_string(index=False))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Update genebuild_status table from production DB.")
    parser.add_argument("-p", "--password", required=True, help="MySQL password for write user")
    parser.add_argument("-t", "--test", action="store_true", help="Run in test mode (no DB updates)")
    parser.add_argument("-am", "--apply_method", action="store_true", help="Apply method updates (default is stop)")
    parser.add_argument("-or", "--old_registry", action="store_true", help="Check old registry status")
    parser.add_argument("-ao", "--apply_old", action="store_true", help="Apply updates from old registry (default is stop).")

    args = parser.parse_args()

    main(args.password, args.test, args.apply_method, args.old_registry, args.apply_old)




