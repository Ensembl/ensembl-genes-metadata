import pymysql
from helper import mysql_fetch_data
import pandas as pd
from logger_settings import get_logger

logger = get_logger(__name__)

def check_status_production_db(gca_tuple):
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