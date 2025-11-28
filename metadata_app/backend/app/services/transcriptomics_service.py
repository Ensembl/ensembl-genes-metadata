# metadata_app/backend/app/services/transcriptomics_service.py

import pandas as pd
import logging
from metadata_app.backend.app.core.database import get_db_connection

def get_transcriptomic_assessment_for_ids(taxon_ids: list) -> pd.DataFrame:
    """
    Retrieve transcriptomic assessment data for a list of taxon IDs.

    Args:
        taxon_ids (list of int): Taxon IDs to retrieve transcriptomic assessment for.

    Returns:
        pd.DataFrame: DataFrame with columns ['taxon_id', 'transc_assess_date', 'aligned_count'].
    """
    if not taxon_ids:
        logging.warning("No taxon IDs provided for transcriptomic assessment.")
        return pd.DataFrame(columns=['taxon_id', 'transc_assess_date', 'aligned_count'])

    try:
        # Remove duplicates and ensure integers
        taxon_ids = list(set(int(tid) for tid in taxon_ids))

        placeholders = ','.join(['%s'] * len(taxon_ids))

        query = f"""
            SELECT m.taxon_id,
                   m.last_check AS transc_assess_date,
                   COUNT(r.qc_status) AS aligned_count
            FROM meta m
            JOIN run r ON m.taxon_id = r.taxon_id
            WHERE m.taxon_id IN ({placeholders})
            GROUP BY m.taxon_id, m.last_check;
        """

        with get_db_connection("transcriptomic") as conn:
            cursor = conn.cursor()
            cursor.execute(query, taxon_ids)
            results = cursor.fetchall()

        df = pd.DataFrame(results, columns=['taxon_id', 'transc_assess_date', 'aligned_count'])
        logging.info(f"Retrieved transcriptomic data for {len(df)} of {len(taxon_ids)} taxon IDs.")

        missing_count = len(taxon_ids) - df['taxon_id'].nunique()
        if missing_count > 0:
            logging.info(f"{missing_count} taxon IDs had no transcriptomic assessment data.")

        return df

    except Exception as e:
        logging.error(f"Error retrieving transcriptomic assessment data: {str(e)}")
        raise


def add_transc_data_to_df(info_df):
    """
    Add transcriptomic assessment data to the information DataFrame.
    """
    try:
        species_ids = info_df['species_taxon_id'].dropna().astype(int).tolist()
        lowest_ids = info_df['lowest_taxon_id'].dropna().astype(int).tolist()
        all_ids = list(set(species_ids + lowest_ids))

        trans_df = get_transcriptomic_assessment_for_ids(all_ids)  # Function accepts list of IDs

        # Merge for species
        info_df = info_df.merge(
            trans_df.rename(columns={
                'transc_assess_date': 'species_transc_assess_date',
                'aligned_count': 'species_aligned_count'
            }),
            how='left',
            left_on='species_taxon_id',
            right_on='taxon_id'
        ).drop(columns=['taxon_id'])

        # Merge for lowest taxon
        info_df = info_df.merge(
            trans_df.rename(columns={
                'transc_assess_date': 'lowest_transc_assess_date',
                'aligned_count': 'lowest_aligned_count'
            }),
            how='left',
            left_on='lowest_taxon_id',
            right_on='taxon_id'
        ).drop(columns=['taxon_id'])

        return info_df

    except Exception as e:
        logging.error(f"Error adding transcriptomic data to DataFrame: {str(e)}")
        raise