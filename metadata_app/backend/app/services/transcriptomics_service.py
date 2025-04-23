# metadata_app/backend/app/services/transcriptomics_service.py

import logging
import pandas as pd
from typing import Dict, List, Set, Optional

from metadata_app.backend.app.core.database import get_db_connection


def get_trancriptomic_assessment(taxonomy_dict):
	"""
	Retrieve transcriptomic assessment data for taxon IDs in the taxonomy dictionary.

	Args:
		taxonomy_dict: Dictionary mapping lowest taxon IDs to their taxonomy hierarchies

	Returns:
		DataFrame with transcriptomic assessment data
	"""
	try:
		# Connect to database
		with get_db_connection("transcriptomic") as conn:
			cursor = conn.cursor()

			logging.info(f"Retrieving trancriptomic assesment data from the registry.")

			trans_taxon_ids = set()

			for lowest_taxon_id, tax_list in taxonomy_dict.items():
				trans_taxon_ids.add(lowest_taxon_id)  # Always include lowest_taxon_id
				for taxon in tax_list:
					if taxon['taxon_class'] in ['species', 'genus']:
						trans_taxon_ids.add(taxon['taxon_class_id'])

			trans_taxon_ids = list(trans_taxon_ids)

			# Build and execute query
			if not trans_taxon_ids:
				logging.warning("No taxon IDs provided for transcriptomic assessment.")
				return pd.DataFrame(columns=['taxon_id', 'transc_assess_date', 'transc_status'])

			trans_placeholders = ','.join(['%s'] * len(trans_taxon_ids))
			where_clause = f"WHERE m.taxon_id IN ({trans_placeholders})"

			query = f"""
                SELECT m.taxon_id, m.last_check AS transc_assess_date, r.qc_status AS transc_status
                FROM meta m
                JOIN run r ON m.taxon_id = r.taxon_id
                {where_clause};
            """

			cursor.execute(query, trans_taxon_ids)
			trans_results = cursor.fetchall()

			# Create DataFrame from results
			trans_df = pd.DataFrame(trans_results, columns=['taxon_id', 'transc_assess_date', 'transc_status'])

			# Log statistics
			found_ids = trans_df['taxon_id'].nunique()
			logging.info(
				f"Transcriptomic assessment data retrieved for {found_ids} taxon IDs out of {len(trans_taxon_ids)} provided.")

			missing_count = len(trans_taxon_ids) - found_ids
			if missing_count > 0:
				logging.info(f"{missing_count} taxon IDs had no transcriptomic assessment data.")

			return trans_df

	except Exception as e:
		logging.error(f"Error retrieving transcriptomic assessment data: {str(e)}")
		raise


def add_transc_data_to_df(info_df, taxonomy_dict):
	"""
	Add transcriptomic assessment data to the information DataFrame.

	Args:
		info_df: DataFrame containing taxonomic information
		taxonomy_dict: Dictionary mapping lowest taxon IDs to their taxonomy hierarchies

	Returns:
		DataFrame with added transcriptomic assessment data
	"""
	try:
		# Get transcriptomic assessment data
		trans_df = get_trancriptomic_assessment(taxonomy_dict)

		# Merge data at each taxonomic level
		for level in ['lowest', 'species', 'genus']:
			merged = info_df[[f"{level}_taxon_id"]].merge(
				trans_df,
				how='left',
				left_on=f"{level}_taxon_id",
				right_on='taxon_id'
			)
			info_df[f"{level}_transc_assess_date"] = merged['transc_assess_date']
			info_df[f"{level}_transc_status"] = merged['transc_status']

		return info_df

	except Exception as e:
		logging.error(f"Error adding transcriptomic data to DataFrame: {str(e)}")
		raise