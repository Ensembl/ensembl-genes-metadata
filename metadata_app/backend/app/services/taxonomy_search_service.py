import logging
from fastapi import HTTPException
from metadata_app.backend.app.core.database import get_db_connection

def search_taxonomy(user_input: str):
    """Search taxonomy by name or ID, return flat list of displayable strings."""
    try:
        with get_db_connection("meta") as conn:
            cursor = conn.cursor()

            query = """
                SELECT DISTINCT t.taxon_class_name, t.taxon_class_id
                FROM taxonomy_name t
                WHERE LOWER(t.taxon_class_name) LIKE LOWER(%s)
                   OR CAST(t.taxon_class_id AS CHAR) LIKE %s
                ORDER BY t.taxon_class_name, t.taxon_class_id
                LIMIT 10;
            """

            like_term = f"%{user_input}%"
            cursor.execute(query, (like_term, like_term))
            rows = cursor.fetchall()
            print("ROWS:", rows)

            if not rows:
                raise HTTPException(
                    status_code=404,
                    detail="No taxonomy found matching your query."
                )

            # Flatten to a set of unique, non-None strings
            results_set = set()
            for row in rows:
                if row["taxon_class_name"]:
                    results_set.add(row["taxon_class_name"])
                if row["taxon_class_id"]:
                    results_set.add(row["taxon_class_id"])
            print(results_set)

            return list(results_set)

    except Exception as e:
        logging.error("Error during taxonomy search", exc_info=True)
        raise HTTPException(
            status_code=500,
            detail=f"Internal server error: {str(e)}"
        )