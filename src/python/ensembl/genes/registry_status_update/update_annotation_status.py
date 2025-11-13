import argparse
import pandas as pd
import pymysql
import logging
from copy_status_old_registry import insert_entries_from_old_registry
from helper import mysql_fetch_data
from production_check import check_status_production_db
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
                gb_status
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

    Args:
        merged_df (pd.DataFrame): Merged DataFrame with columns:
            - 'genebuild_status_id'
            - 'gb_status' (current registry status)
            - 'status' (production status)
            - 'release_date' (current registry release date)
            - 'release_date_production' (from production)
            - 'genebuild_version'

    Returns:
        pd.DataFrame: DataFrame with updates where status or release date changes.
    """
    df = merged_df.copy()

    # Normalize string values
    df['gb_status'] = df['gb_status'].astype(str).str.strip().str.lower()
    df['status'] = df['status'].astype(str).str.strip().str.capitalize()

    # Normalize release dates (convert to datetime)
    df['release_date'] = pd.to_datetime(df.get('release_date', pd.NaT), errors='coerce')
    df['release_date_production'] = pd.to_datetime(df.get('release_date_production', df.get('release_date')), errors='coerce')

    # Initialize new status and release date columns
    df['gb_status_new'] = df['gb_status']
    df['release_date_new'] = df['release_date']

    # 1. Skip faulty entries
    condition_faulty = df['status'] == 'Faulty'
    if condition_faulty.any():
        logger.info(f"Skipping {condition_faulty.sum()} 'Faulty' entries.")

        # Keep only rows with faulty status
        faulty_ds = df[condition_faulty].copy()

        # Save faulty GCAs to file
        faulty_ds.to_csv("faulty_status.csv", index=False)
        logger.info("Saved faulty production status to 'faulty_status.csv'.")

        # Remove rows with faulty status from the original df
        df = df[~condition_faulty].copy()
        logger.info(f"{len(df)} rows remain after removing entries with faulty status.")


    # 2. Released in production but not 'live' in registry
    condition_released = (df['status'] == 'Released') & (df['gb_status'] != 'live')
    df.loc[condition_released, 'gb_status_new'] = 'live'
    df.loc[condition_released, 'release_date_new'] = df.loc[condition_released, 'release_date_production']

    # 3. If already 'live' but missing release date → fill it from production
    condition_missing_date = (df['gb_status'] == 'live') & (df['release_date'].isna()) & (df['release_date_production'].notna())
    df.loc[condition_missing_date, 'release_date_new'] = df.loc[condition_missing_date, 'release_date_production']

    # 4. If already 'live' but release date mismatch → update it
    condition_mismatch_date = (
        (df['gb_status'] == 'live')
        & (df['release_date'].notna())
        & (df['release_date_production'].notna())
        & (df['release_date'] != df['release_date_production'])
    )
    df.loc[condition_mismatch_date, 'release_date_new'] = df.loc[condition_mismatch_date, 'release_date_production']

    # 5. Processing / Submitted → handed_over
    condition_handed_over = df['status'].isin(['Processed', 'Processing', 'Submitted']) & (
        df['gb_status'] != 'handed_over'
    )
    df.loc[condition_handed_over, 'gb_status_new'] = 'handed_over'

    # Compare status
    status_changed = df['gb_status'] != df['gb_status_new']

    # Compare release date, treating NaT properly
    release_changed = ~((df['release_date_new'].isna() & df['release_date'].isna()) |
                        (df['release_date_new'] == df['release_date']))

    # Keep only rows where either changed
    updated_df = df[status_changed | release_changed].copy()

    updated_df = updated_df[['genebuild_status_id', 'gb_status_new', 'release_date_new', 'genebuild_version']]
    logger.info(f"Found {len(updated_df)} annotations requiring status or release date updates.")

    return updated_df


def update_genebuild_status(updated_df, password):
    """
    Update the genebuild_status table with new statuses and release dates.

    Args:
        updated_df (pd.DataFrame): DataFrame with columns:
            - genebuild_status_id
            - gb_status_new
            - release_date_new (optional, only for 'live')
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

    try:
        with connection.cursor() as cursor:
            for _, row in updated_df.iterrows():
                status = row['gb_status_new']
                genebuild_status_id = row['genebuild_status_id']
                release_date = row.get('release_date_new', None)

                if status == 'live' and pd.notnull(release_date):
                    sql = """
                        UPDATE genebuild_status
                        SET gb_status = %s,
                            release_date = %s
                        WHERE genebuild_status_id = %s
                    """
                    cursor.execute(sql, (status, release_date, genebuild_status_id))
                else:
                    sql = """
                        UPDATE genebuild_status
                        SET gb_status = %s
                        WHERE genebuild_status_id = %s
                    """
                    cursor.execute(sql, (status, genebuild_status_id))

        logging.info(f"Updated {len(updated_df)} genebuild_status rows.")

    except pymysql.Error as e:
        logging.error("MySQL error: %s", e)
        raise

    finally:
        connection.close()


def main(password, test, old_registry, stop_appy):

    if old_registry:
        logger.info("Checking status in old registry.")
        gb_status = get_genebuild_status()
        copy_from_old_registry = insert_entries_from_old_registry(password, gb_status, stop_appy)


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
    production_status = check_status_production_db(gca_tuple)

    # Ensure consistent column names
    if 'gca_accession' not in production_status.columns:
        logger.error("Expected columns missing from production DB query.")
        return

    # Merge the two dataframes on 'gca_accession'
    merged_df = pd.merge(
        gb_status,
        production_status,
        on='gca_accession',
        how='inner',
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
    parser.add_argument("-or", "--old_registry", type=bool, help="Check old registry status.")
    parser.add_argument("-sa", "--stop_apply", type=bool, default=True, help="If true don't apply updates from old registry.")


    args = parser.parse_args()

    main(args.password, args.test, args.old_registry, args.stop_apply)




