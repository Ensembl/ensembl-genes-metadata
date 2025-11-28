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
                    CONCAT(a.gca_chain, ".", a.gca_version) AS gca
                FROM genebuild_status g
                JOIN assembly a ON a.assembly_id = g.assembly_id
                LEFT JOIN bioproject b ON b.assembly_id = a.assembly_id
                LEFT JOIN main_bioproject m ON m.bioproject_id = b.bioproject_id
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

        collapsed = (
            df.groupby("gca", as_index=False)
            .agg({
                # join all bioproject names as a single string
                "bioproject_name": lambda x: ", ".join(sorted(set(x))),
                # for all other columns, take any value (they are identical)
                # replace these with your actual column names
                "genebuild_status_id": "first",
                "date_status_update": "first",
                "genebuilder": "first",
                "gb_status": "first",
            })
        )




        current_genebuilders = ["lazar", "leanne", "jackt", "vianey", "swati", "ereboperezsilva", "ftricomi"]

        df_filtered = collapsed[collapsed["genebuilder"].isin(current_genebuilders)]

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