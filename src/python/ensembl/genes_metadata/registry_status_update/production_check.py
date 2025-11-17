"""
production_check.py

This module provides functionality to fetch the genebuild status of genome assemblies
from the Ensembl production database. It queries multiple tables to retrieve the
current genebuild status, release date, and genebuild version for a list of assemblies.

Functions:
----------
check_status_production_db(gca_tuple)
    Retrieves production genebuild status information for the specified GCAs.

Dependencies:
-------------
- pandas
- pymysql
- helper (module providing mysql_fetch_data function)
- logger_settings (module providing get_logger function)
"""

import pymysql
from helper import mysql_fetch_data
import pandas as pd
from logger_settings import get_logger

logger = get_logger(__name__)

def check_status_production_db(gca_tuple):
    """
    Fetch genebuild status, release date, and genebuild version for given assemblies
    from the Ensembl production database.

    Args:
        gca_tuple (tuple): Tuple of GCA accession strings to query.

    Returns:
        pd.DataFrame: DataFrame containing the following columns:
            - gca_accession: Assembly accession (GCA).
            - status: Genebuild status in production (e.g., Released, Processed).
            - release_date: Ensembl release date for the genebuild.
            - genebuild_version: Version of the genebuild dataset.

    Notes:
        - Only includes datasets with dataset.name = "genebuild".
        - Only considers current genome datasets (genome_dataset.is_current = 1).
        - Drops duplicate GCA entries, keeping the first occurrence.
        - Logs the number of entries found.
        - Returns an empty DataFrame in case of MySQL errors.
    """
    try:
        production_query = f"""
            SELECT 
                assembly.accession AS gca_accession,
                dataset.status,
                ensembl_release.release_date,
                dataset_attribute.value AS genebuild_version
            FROM assembly
            LEFT JOIN genome
                ON assembly.assembly_id = genome.assembly_id
            LEFT JOIN genome_dataset
                ON genome.genome_id = genome_dataset.genome_id
            LEFT JOIN dataset
                ON genome_dataset.dataset_id = dataset.dataset_id
            LEFT JOIN dataset_attribute
                ON dataset.dataset_id = dataset_attribute.dataset_id
            LEFT JOIN dataset_source
                ON dataset.dataset_source_id = dataset_source.dataset_source_id
            LEFT JOIN ensembl_release on genome_dataset.release_id = ensembl_release.release_id
            WHERE dataset.name = "genebuild"
                AND assembly.accession IN {gca_tuple}
                AND genome_dataset.is_current = 1
                AND dataset_attribute.attribute_id = 71
        """

        production_status = mysql_fetch_data(
            production_query,
            host="mysql-ens-production-1",
            user="ensro",
            port=4721,
            database="ensembl_genome_metadata",
            password=""
        )
        production_status = pd.DataFrame(production_status)
        production_status = production_status.drop_duplicates(subset='gca_accession', keep='first')

        logger.info(f"Found {len(production_status)} entries in production table.")
        return production_status

    except pymysql.Error as err:
        logger.error("MySQL error: %s", err)
        return