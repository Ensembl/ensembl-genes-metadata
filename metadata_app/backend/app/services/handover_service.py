# app/services/home_page_service.py
import logging

import numpy as np
import pandas as pd
from metadata_app.backend.app.core.database import get_db_connection
import pymysql.cursors


def get_ready_to_ho():
    try:
        with get_db_connection("meta") as conn:
            cursor = conn.cursor(pymysql.cursors.DictCursor)  # fetch rows as dicts
            query = """
                SELECT 
                    g.genebuild_status_id,
                    g.gb_status,
                    g.date_status_update,
                    g.genebuilder, 
                    m.bioproject_name,
                    s.scientific_name,
                    CONCAT(a.gca_chain, ".", a.gca_version) AS gca
                FROM genebuild_status g
                JOIN assembly a ON a.assembly_id = g.assembly_id
                LEFT JOIN bioproject b ON b.assembly_id = a.assembly_id
                LEFT JOIN main_bioproject m ON m.bioproject_id = b.bioproject_id
                LEFT JOIN species s ON s.lowest_taxon_id = a.lowest_taxon_id
                WHERE g.gb_status IN ('completed', 'pre_released')
            """
            cursor.execute(query)
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
                "genebuilder": "first",
                "gb_status": "first",
                "scientific_name": "first",
            })
        )

        current_genebuilders = ["lazar", "leanne", "jackt", "vianey", "swati", "ereboperezsilva", "ftricomi"]

        df_filtered = collapsed[collapsed["genebuilder"].isin(current_genebuilders)]

        df_filtered["production_name"] = (
                df_filtered["scientific_name"]
                .str.lower()
                .str.replace(" ", "_")
                + "_"
                + df_filtered["gca"]
                .str.lower()
                .str.replace(".", "v", regex=False)
                .str.replace("_", "", regex=False)
        )

        rename_map = {
            "lazar": "Anna",
            "leanne": "Leanne",
            "jackt": "Jack",
            "vianey": "Vianey",
            "swati": "Swati",
            "ereboperezsilva": "Jose",
            "ftricomi": "Francesca",
        }

        grouped = {
            rename_map.get(genebuilder, genebuilder): rows.drop(columns=['genebuilder']).to_dict(orient="records")
            for genebuilder, rows in df_filtered.groupby("genebuilder")
        }

        return grouped

    except Exception as e:
        logging.error(f"Error in get_ready_to_ho: {e}", exc_info=True)
        return {}