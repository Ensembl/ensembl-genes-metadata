# app/services/home_page_service.py
import logging

import numpy as np
import pandas as pd
from metadata_app.backend.app.core.database import get_db_connection
import pymysql.cursors


def get_ready_to_ho(genebuilder):
	try:
		with get_db_connection("meta") as conn:
			cursor = conn.cursor(pymysql.cursors.DictCursor)

			query = """
				SELECT 
					g.genebuild_status_id,
					g.gb_status,
					g.date_status_update,
					m.bioproject_name,
					s.scientific_name,
					g.annotation_method,
					CONCAT(a.gca_chain, ".", a.gca_version) AS gca
				FROM genebuild_status g
				JOIN assembly a ON a.assembly_id = g.assembly_id
				LEFT JOIN bioproject b ON b.assembly_id = a.assembly_id
				LEFT JOIN main_bioproject m ON m.bioproject_id = b.bioproject_id
				LEFT JOIN species s ON s.lowest_taxon_id = a.lowest_taxon_id
				WHERE g.gb_status IN ('completed', 'pre_released', 'handed_over', 'in_progress', 'check_busco', 'insufficient_data')
				  AND g.genebuilder = %s
			"""

			cursor.execute(query, (genebuilder,))
			result = cursor.fetchall()

		if not result:
			return {}


		# Replace NaN/Inf with the string "None" so JSON serialization is safe
		# After building df
		df = pd.DataFrame(result)

		# Replace NaN/inf with Unknown
		df = df.replace([np.nan, np.inf, -np.inf], "Unknown")

		def clean_bioproject_names(values):
			unique_vals = {v for v in values if v != "Unknown"}
			if unique_vals:
				# real values exist → return only the real ones
				return ", ".join(sorted(unique_vals))
			else:
				# only "Unknown" present
				return "Unknown"

		collapsed = (
			df.groupby("gca", as_index=False)
			.agg({
				"bioproject_name": clean_bioproject_names,
				"genebuild_status_id": "first",
				"date_status_update": "first",
				"gb_status": "first",
				"scientific_name": "first",
				"annotation_method": "first",
			})
		)

		ho_ready = ["pre_released", "completed"]

		df_filtered = collapsed[collapsed["gb_status"].isin(ho_ready)]

		df_ready = df_filtered

		df_ready.loc[:, "production_name"] = (
				df_ready["scientific_name"]
				.str.lower()
				.str.replace(" ", "_")
				+ "_"
				+ df_ready["gca"]
				.str.lower()
				.str.replace(".", "v", regex=False)
				.str.replace("_", "", regex=False)
		)

		# check data and in progress for more than 6 months
		date_series = pd.to_datetime(collapsed["date_status_update"], errors="coerce")
		six_months_ago = pd.Timestamp.today() - pd.DateOffset(months=6)

		count_ho_ready = df_filtered["gca"].nunique()

		# handed over (from collapsed)
		check_data = ["check_busco", "insufficient_data"]
		check_data_df = collapsed[collapsed["gb_status"].isin(check_data)]

		count_data = (
			check_data_df.loc[
				(date_series < six_months_ago),
				"gca"
			]
		).nunique()

		list_data = (
			check_data_df.loc[
				(date_series < six_months_ago),
				["gca", "scientific_name", "gb_status"]
			]
		)

		# pending
		count_pending = (
			collapsed.loc[
				(collapsed["gb_status"] == "in_progress") &
				(date_series < six_months_ago),
				"gca"
			]
		).nunique()

		list_pending = (
			collapsed.loc[
				(collapsed["gb_status"] == "in_progress") &
				(date_series < six_months_ago),
				["gca", "scientific_name", "gb_status"]
			]
		)
		print("Count pending:", count_pending)
		print("List data:", list_data)

		return (
			df_ready.to_dict(orient="records"),
			count_ho_ready,
			count_data,
			count_pending,
			pd.DataFrame(list_data).to_dict(orient="records") if not list_data.empty else [],
			pd.DataFrame(list_pending).to_dict(orient="records") if not list_pending.empty else []
		)

	except Exception as e:
		logging.error(f"Error in get_ready_to_ho: {e}", exc_info=True)
		return {}

def update_gca(genebuilder, items):
	try:
		if not items:
			return 0

		with get_db_connection("meta_write") as conn:
			cursor = conn.cursor(pymysql.cursors.DictCursor)

			query = """
				UPDATE
					genebuild_status g
				SET g.gb_status = %s,
                    g.date_status_update = CURRENT_DATE
				WHERE g.gb_status IN ('completed', 'pre_released')
				  AND g.genebuilder = %s
				  AND g.gca_accession = %s
				  AND g.annotation_method = %s
			"""
			updated_rows = 0

			for item in items:
				cursor.execute(
					query,
					(
						"abandoned",
						genebuilder,
						item["gca"],
						item["annotation_method"],
					),
				)
				updated_rows += cursor.rowcount
				if cursor.rowcount == 0:
					logging.warning(
						f"No rows updated for GCA={item['gca']} "
						f"method={item['annotation_method']} "
						f"genebuilder={genebuilder}"
					)

			conn.commit()
			return updated_rows


	except Exception as e:
		logging.error(f"Error in update_gca: {e}", exc_info=True)
		return 0
