# app/services/home_page_service.py

import json
import logging
from metadata_app.backend.app.core.database import get_db_connection
import pandas as pd
import datetime
from metadata_app.backend.app.services.taxonomy_service import get_descendant_taxa


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

def get_annotations_per_genebuilder():
    try:
        with get_db_connection("meta") as conn:
            cursor = conn.cursor()
            query = """
                SELECT g.genebuild_id, g.date_completed, g.genebuilder
                FROM genebuild g
            """
            cursor.execute(query)
            result = cursor.fetchall()

        df = pd.DataFrame(result, columns=["genebuild_id", "date_completed", "genebuilder"])

        if df.empty:
            return []

        # Filter to last 30 days
        logging.info(f"getting timeframe")
        thirty_days_ago = (datetime.now() - datetime.timedelta(days=30)).date()
        df["date_completed"] = pd.to_datetime(df["date_completed"]).dt.date
        recent_df = df[df["date_completed"] >= thirty_days_ago]
        logging.info(f"filtered last 30 days")
        logging.debug(print(recent_df))

        # Group by genebuilder and count
        count_df = recent_df.groupby("genebuilder").size().reset_index(name="annotations")

        return count_df.to_dict(orient='records')

    except Exception as e:
        # Optionally log or re-raise
        print(f"Error in get_annotations_per_genebuilder: {e}")
        return []



def bin_by_genebuild_method(bioproject_id, release_type, taxon_id, release_date):
    """Bins assemblies based on genebuild.method."""
    try:
        with get_db_connection("meta") as conn:
            cursor = conn.cursor()
            # Build query conditions
            conditions = []
            params = []

            if bioproject_id:
                conditions.append(f"b.bioproject_id IN ({','.join(['%s'] * len(bioproject_id))})")
                params.extend(bioproject_id)
                logging.info(f"Filtering by BioProject IDs: {', '.join(bioproject_id)}")

            if release_date:
                if isinstance(release_date, pd.Timestamp):
                    release_date = release_date.strftime('%Y-%m-%d')
                elif isinstance(release_date, (datetime.date, datetime.datetime)):
                    release_date = release_date.strftime('%Y-%m-%d')
                conditions.append("g.release_date_beta >= %s")
                params.append(release_date)
                logging.info(f"Filtering by release date: {release_date}")

            if taxon_id:
                descendant_taxa = get_descendant_taxa(taxon_id)
                if not descendant_taxa:
                    logging.error(f"No descendant taxon IDs found for {taxon_id}.")
                    return f"No descendant taxa found for Taxon ID {taxon_id}.", None, None, None, None
                conditions.append(f"a.lowest_taxon_id IN ({','.join(['%s'] * len(descendant_taxa))})")
                params.extend(descendant_taxa)
                logging.info(f"Filtering by lowest taxon ID: {', '.join(str(id) for id in descendant_taxa)}")

            if release_type:
                conditions.append(f"g.release_type IN ({','.join(['%s'] * len(release_type))})")
                params.extend(release_type)
                logging.info(f"Filtering by Release Type: {', '.join(release_type)}")

            # If there are conditions, join them with AND; otherwise, select all
            where_clause = " WHERE " + " AND ".join(conditions) if conditions else ""

            query = f"""
                    SELECT g.annotation_method, g.release_type, g.release_date_beta, b.bioproject_id, a.lowest_taxon_id
                    FROM genebuild g
                    LEFT JOIN assembly a ON g.assembly_id = a.assembly_id
                    LEFT JOIN bioproject b on g.assembly_id = b.assembly_id
                    {where_clause};
                """
            cursor.execute(query, params)
            result = cursor.fetchall()

        df = pd.DataFrame(result, columns=["annotation_method", "release_type", "release_date_beta", "bioproject_id", "lowest_taxon_id"])

        if df.empty:
            return []


        # Group by genebuilder and count
        method_summary = df.groupby('annotation_method', observed=False).size().reset_index(name='number_of_annotations')

        return method_summary.to_dict(orient='records')

    except Exception as e:
        # Optionally log or re-raise
        print(f"Error in bin_by_genebuld_method: {e}")
        return []

