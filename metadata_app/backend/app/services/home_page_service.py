# app/services/home_page_service.py

import json
import logging
from metadata_app.backend.app.core.database import get_db_connection
import pandas as pd


def load_bioproject_mapping():
    """Hardcoded path for clade settings."""
    json_file = "data/bioproject_mapping.json"
    with open(json_file, "r") as f:
        logging.info("Loading bioproject mapping json file.")
        return json.load(f)


def get_annotation_counts_by_bioproject():
    """Returns a count of annotations per mapped BioProject, with friendly names."""
    try:
        with get_db_connection("meta") as conn:
            cursor = conn.cursor()
            query = """
                SELECT b.bioproject_id, COUNT(g.assembly_id) AS annotation_count
				FROM bioproject b
				JOIN genebuild g ON b.assembly_id = g.assembly_id
				WHERE g.gb_status != 'in_progress'
				GROUP BY b.bioproject_id;
            """
            cursor.execute(query)
            result = cursor.fetchall()

        df = pd.DataFrame(result)
        if df.empty:
            return []

        mapping = load_bioproject_mapping()

        # Only keep mapped BioProjects
        df = df[df["bioproject_id"].isin(mapping.keys())]
        df["associated_project"] = df["bioproject_id"].map(mapping)

        return df[["bioproject_id", "associated_project", "annotation_count"]].to_dict(orient="records")

    except Exception as e:
        logging.error(f"Error fetching annotation counts: {e}")
        return []


def get_assemblies_per_year():
	try:
		with get_db_connection("meta") as conn:
			cursor = conn.cursor()
			query = """
                SELECT a.assembly_id, a.release_date
                FROM assembly a
                WHERE a.is_current = 'current'
            """
			cursor.execute(query)
			result = cursor.fetchall()

		df = pd.DataFrame(result, columns=["assembly_id", "release_date"])
		if df.empty:
			return []

		# Convert release_date to datetime and extract year
		df['release_date'] = pd.to_datetime(df['release_date'], errors='coerce')
		df['year'] = df['release_date'].dt.year

		# Filter to include only 2019 and later
		df = df[df['year'] >= 2019]

		# Group by year and count
		counts_by_year = df.groupby('year')['assembly_id'].count().reset_index()
		counts_by_year.columns = ['year', 'assembly_count']

		return counts_by_year.to_dict(orient='records')

	except Exception as e:
		print(f"Error in get_assemblies_per_year: {e}")
		return []

def get_annotations_per_year():
	try:
		with get_db_connection("meta") as conn:
			cursor = conn.cursor()
			query = """
                SELECT g.genebuild_id, g.date_completed
                FROM genebuild g
            """
			cursor.execute(query)
			result = cursor.fetchall()

		df = pd.DataFrame(result, columns=["genebuild_id", "date_completed"])
		if df.empty:
			return []

		# Convert date_completed to datetime and extract year
		df['date_completed'] = pd.to_datetime(df['date_completed'], errors='coerce')
		df['year'] = df['date_completed'].dt.year
		# Filter to include only 2019 and later
		df = df[df['year'] >= 2019]

		# Group by year and count
		counts_by_year = df.groupby('year')['genebuild_id'].count().reset_index()
		counts_by_year.columns = ['year', 'annotation_count']

		return counts_by_year.to_dict(orient='records')

	except Exception as e:
		print(f"Error in get_annotations_per_year: {e}")
		return []

def get_metadata_registry_update_dates():
	try:
		with get_db_connection("meta") as conn:
			cursor = conn.cursor()
			query = """
                SELECT date_value
                FROM update_date
                WHERE update_type = 'regular_update'
            """
			cursor.execute(query)
			result = cursor.fetchall()

		# Convert result to DataFrame with correct column name
		df = pd.DataFrame(result, columns=["date_value"])

		if df.empty:
			return []

		# Convert to string list of dates
		return df["date_value"].astype(str).tolist()

	except Exception as e:
		print(f"Error fetching metadata registry update dates: {e}")
		return []

def get_transcriptomic_registry_update_dates():
	try:
		with get_db_connection("transcriptomic") as conn:
			cursor = conn.cursor()
			query = """
                SELECT last_check
                FROM meta
            """
			cursor.execute(query)
			result = cursor.fetchall()

		# Convert result to DataFrame with correct column name
		df = pd.DataFrame(result, columns=["last_check"])

		if df.empty:
			return []

		# Convert to string list of dates
		return df["last_check"].astype(str).tolist()

	except Exception as e:
		print(f"Error fetching transcriptomic registry update dates: {e}")
		return []
