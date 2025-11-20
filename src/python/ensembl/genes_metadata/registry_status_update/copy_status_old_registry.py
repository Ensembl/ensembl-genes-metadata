"""
copy_status_old_registry.py

This module provides functionality to migrate and synchronize genebuild status entries
from the old genebuild registry to the new registry. It handles fetching, filtering,
mapping, and inserting old entries, ensuring only valid and current datasets are applied.

Key Features:
--------------
- Fetch the latest genebuild status entries from the old registry.
- Filter out entries that do not exist in the production database.
- Map old registry statuses to the new registry format.
- Retrieve corresponding assembly IDs from the new registry.
- Insert missing entries into the new genebuild_status table.
- Supports a test mode to preview changes without applying them.

Dependencies:
-------------
- numpy
- pandas
- pymysql
- helper (module providing mysql_fetch_data)
- production_check (module providing check_status_production_db)
- logger_settings (module providing get_logger)
"""
import numpy as np
import pandas as pd
import pymysql
from helper import mysql_fetch_data
from production_check import check_status_production_db
from logger_settings import get_logger

logger = get_logger(__name__)



def fetch_status_old_registry():
    """Fetch latest genebuild_status entries from the registry."""

    query = f"""
    SELECT
        genebuild_status.assembly_accession AS gca_accession, 
        genebuild_status.progress_status AS gb_status, 
        genebuild_status.date_started, 
        genebuild_status.date_completed, 
        genebuild_status.genebuilder, 
        genebuild_status.annotation_source
    FROM genebuild_status
    WHERE genebuild_status.is_current = 1
    """

    result = mysql_fetch_data(
        query,
        host="mysql-ens-genebuild-prod-1",
        user="ensro",
        port=4527,
        database="gb_assembly_registry",
        password=""
    )

    old_gb_status = pd.DataFrame(result)

    old_gb_status['date_started'] = pd.to_datetime(old_gb_status['date_started'])
    logger.info(f"Found {len(old_gb_status)} entries in the old genebuild_status table.")
    old_reg_df_latest = (
        old_gb_status
        .sort_values('date_started', ascending=False)
        .drop_duplicates(subset='gca_accession', keep='first')
        .reset_index(drop=True)
    )
    logger.info(f"Keeping only the latest entries from old registry per GCA. Remaining entries: {len(old_reg_df_latest)}")

    return old_reg_df_latest

def delete_entries_not_production_db(filtered_old):
    """Keep all entries except those marked 'live' or 'handed_over' that are NOT found in production db."""
    gca_tuple = tuple(filtered_old['gca_accession'].unique())
    if len(gca_tuple) == 1:
        gca_tuple = (gca_tuple[0],)

    logger.info("Getting GCA status from production DB.")
    production_found = check_status_production_db(gca_tuple)

    # Accessions present in production
    if production_found is None or len(production_found) == 0:
        logger.warning("No accessions found in production DB. Treating all as missing.")
        existing_accessions_prod = set()
    else:
        existing_accessions_prod = set(production_found['gca_accession'])

    # Log counts
    logger.info(f"Before filtering: {len(filtered_old)} rows in old registry")

    # Remove rows that are 'live' or 'handed_over' and not found in production
    not_in_prod = (
        (filtered_old["gb_status"].isin(["live", "handed_over"]))
        & (~filtered_old["gca_accession"].isin(existing_accessions_prod))
    )
    filtered_final = filtered_old[~not_in_prod].copy()

    logger.info(
        f"Removed {not_in_prod.sum()} 'live' or 'handed_over' entries not found in production. "
        f"Final count: {len(filtered_final)} rows."
    )

    return filtered_final


def check_status_update_old_registry(gb_status, old_reg_df_latest):
    """Detect cases where the old registry progressed but the new registry is still in progress."""

    status_map = {
        "handed over": "live",
        "completed": "completed",
        "Completed": "completed",
        "Check BUSCO": "check_busco",
        "BUSCO Check": "check_busco",
        "in progress": "in_progress",
        "Pre-Released": "pre_released",
        "Insufficient Data": "insufficient_data"
    }

    # New registry "in-progress" statuses
    in_progress_set = {
        "in_progress",
        "pre_released",
        "insufficient_data",
    }

    # Old registry "progressed" statuses
    progressed_set = {
        "completed",
        "check_busco",
        "insufficient_data",
        "pre_released",
    }

    # --- 1. Map old statuses ---
    old_reg_df_latest["mapped_status"] = (
        old_reg_df_latest["gb_status"]
        .map(status_map)
        .fillna(old_reg_df_latest["gb_status"])
    )

    # --- 2. Merge with new registry ---
    merged = old_reg_df_latest.merge(
        gb_status[["gca_accession", "gb_status"]],
        on="gca_accession",
        how="left",
        suffixes=("_old", "_new")
    )

    # Map new registry statuses too
    merged["gb_status_new_mapped"] = (
        merged["gb_status_new"]
        .map(status_map)
        .fillna(merged["gb_status_new"])
    )

    # --- 3. Apply final filtering ---
    filtered = merged[
        merged["mapped_status"].isin(progressed_set) &  # old registry progressed
        merged["gb_status_new_mapped"].isin(in_progress_set) &  # new registry still in progress
        (merged["mapped_status"] != merged["gb_status_new_mapped"])  # statuses must differ
        ].copy()

    # --- Logging ---
    for _, row in filtered.iterrows():
        logger.info(
            f"Old registry progressed but new registry did not for {row['gca_accession']}: "
            f"old={row['mapped_status']} new={row['gb_status_new_mapped']}"
        )

    return filtered


def add_entries_from_old_registry(gb_status):
    """Add annotation statuses from old registry to df (filtered and mapped)."""
    old_reg_df_latest = fetch_status_old_registry()
    old_reg_df_latest['gca_accession'] = old_reg_df_latest['gca_accession'].str.strip()
    gb_status['gca_accession'] = gb_status['gca_accession'].str.strip()

    # Check if old registry status changed since copy
    changed = check_status_update_old_registry(gb_status, old_reg_df_latest)

    # Accessions already in gb_status
    existing_accessions = set(gb_status['gca_accession'])

    # Keep only rows not in gb_status
    logger.info(f"Found {len(old_reg_df_latest)} in the old registry")
    logger.info(f"Found {len(existing_accessions)} entries in the new genebuild_status table.")
    filtered_old = old_reg_df_latest[~old_reg_df_latest['gca_accession'].isin(existing_accessions)].copy()
    logger.info(f"Filtered old registry: {len(filtered_old)} rows remain after removing existing accessions")

    # Transform statuses to match new registry
    status_map = {
        "handed over": "live",
        "completed": "completed",
	    "Completed": "completed",
        "Check BUSCO": "check_busco",
        "BUSCO Check": "check_busco",
        "in progress": "in_progress",
        "Pre-Released": "pre_released",
        "Insufficient Data": "insufficient_data"
    }

    # Map the progress_status column to new gb_status
    filtered_old['gb_status'] = filtered_old['gb_status'].map(status_map).fillna(filtered_old['gb_status'])

    # Rename columns
    filtered_old = filtered_old.rename(
        columns={"annotation_source": "annotation_method", "date_completed": "date_status_update"})
    filtered_old["last_genebuild_update"] = filtered_old["date_status_update"]

    filtered_old = delete_entries_not_production_db(filtered_old)

    return filtered_old


def find_assembly_id(gb_status):
    filtered_old = add_entries_from_old_registry(gb_status)

    logger.info("Finding assembly id in new registry.")
    gca_df = filtered_old[['gca_accession']].copy()

    # Convert GCAs to tuple for SQL IN clause
    gca_tuple = tuple(gca_df['gca_accession'].unique())
    if len(gca_tuple) == 1:
        # Single value needs a trailing comma
        gca_tuple = (gca_tuple[0],)

    # SQL query to fetch assembly_id for these GCAs
    query = f"""
    SELECT
        assembly.assembly_id,
            CONCAT(assembly.gca_chain, '.', assembly.gca_version) AS gca_accession
    FROM assembly
    WHERE CONCAT(assembly.gca_chain, '.', assembly.gca_version) IN {gca_tuple}
    """

    result_ai = mysql_fetch_data(
        query,
        host="mysql-ens-genebuild-prod-1",
        user="ensro",
        port=4527,
        database="gb_assembly_metadata",
        password=""
    )

    assembly_id = pd.DataFrame(result_ai)
    logger.info(f"Found {len(assembly_id)} assembly ids found in new the registry.")

    # Merge back to df
    # Left merge with df to keep all rows
    assembly_id_merged = filtered_old.merge(
        assembly_id,
        on='gca_accession',  # both DataFrames have the same column name now
        how='left',
        suffixes=('', '_assembly')  # keeps deduplicated_df columns intact
    )

    assembly_id_merged["last_attempt"] = 1
    assembly_id_merged["annotation_source"] = "ensembl"

    new_columns = ["annotation_method", "genebuild_version", "release_date", "release_type"]

    for col in new_columns:
        assembly_id_merged[col] = None

    return assembly_id_merged



def insert_entries_from_old_registry(password, gb_status, stop_apply):
    assembly_id_merged = find_assembly_id(gb_status)

    # Keep only rows with missing assembly_id
    missing_assembly_df = assembly_id_merged[assembly_id_merged['assembly_id'].isnull()][['gca_accession']].copy()
    logger.info(f"Found {len(missing_assembly_df)} entries with missing assembly_id.")

    # Save missing GCAs to file
    missing_assembly_df.to_csv("missing_assembly_ids.csv", index=False)
    logger.info("Saved missing GCA accessions to 'missing_assembly_ids.csv'.")

    # Remove rows with missing assembly_id from the original df
    assembly_id_merged = assembly_id_merged[assembly_id_merged['assembly_id'].notnull()].copy()
    logger.info(f"{len(assembly_id_merged)} rows remain after removing entries with missing assembly_id.")

    if stop_apply:
        logger.info("Test mode: no updates applied from old registry.")
        logger.info("Annotations that would be inserted to new registry from old:\n%s", assembly_id_merged.to_string(index=False))
        return assembly_id_merged

    # Replace NaN with None (so MySQL will accept them as NULL)
    insert_merged_df = assembly_id_merged.replace({np.nan: None})

    # Build INSERT query
    columns = insert_merged_df.columns.tolist()
    placeholders = ", ".join(["%s"] * len(columns))
    columns_str = ", ".join(columns)

    insert_query = f"INSERT INTO genebuild_status ({columns_str}) VALUES ({placeholders})"

    # Convert DataFrame to list of tuples
    data_to_insert = [tuple(x) for x in insert_merged_df.to_numpy()]

    # Insert into MySQL
    connection = pymysql.connect(host="mysql-ens-genebuild-prod-1", user="ensadmin", password=password, port=4527, database="gb_assembly_metadata")
    try:
        with connection.cursor() as cursor:
            cursor.executemany(insert_query, data_to_insert)
        connection.commit()
    finally:
        connection.close()

    logger.info(f"Inserted {len(insert_merged_df)} rows from old registry to new registry.")

    return assembly_id_merged, update_from_old_registry

