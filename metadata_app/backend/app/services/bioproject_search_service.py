import logging
from fastapi import HTTPException
from metadata_app.backend.app.core.database import get_db_connection

def search_bioproject(user_input: str):
    """Search BioProject by ID or name (fuzzy match), return flat list of displayable strings."""
    try:
        with get_db_connection("meta") as conn:
            cursor = conn.cursor()

            query = """
                SELECT DISTINCT b.bioproject_id, mb.bioproject_name
                FROM bioproject b
                LEFT JOIN main_bioproject mb ON mb.bioproject_id = b.bioproject_id
                WHERE LOWER(mb.bioproject_name) LIKE LOWER(%s)
                   OR LOWER(b.bioproject_id) LIKE LOWER(%s)
                ORDER BY mb.bioproject_name, b.bioproject_id
                LIMIT 10;
            """

            like_term = f"%{user_input}%"
            cursor.execute(query, (like_term, like_term))
            rows = cursor.fetchall()
            print("ROWS:", rows)

            if not rows:
                raise HTTPException(
                    status_code=404,
                    detail="No BioProjects found matching your query."
                )

            # Flatten to a set of unique, non-None strings
            results_set = set()
            for row in rows:
                if row["bioproject_id"]:
                    results_set.add(row["bioproject_id"])
                if row["bioproject_name"]:
                    results_set.add(row["bioproject_name"])
            print(results_set)

            return list(results_set)

    except Exception as e:
        logging.error("Error during BioProject search", exc_info=True)
        raise HTTPException(
            status_code=500,
            detail=f"Internal server error: {str(e)}"
        )