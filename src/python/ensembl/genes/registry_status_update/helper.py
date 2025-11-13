import pymysql
from logger_settings import get_logger

logger = get_logger(__name__)

def mysql_fetch_data(query, host, user, password, port, database):
    """
    Execute a MySQL query and return the results as a list of dictionaries.

    Args:
        query (str): SQL query to execute.
        host (str): MySQL host.
        user (str): MySQL username.
        password (str): MySQL password.
        port (int): MySQL port.
        database (str): Database name.

    Returns:
        list[dict]: List of rows as dictionaries.
    """
    try:
        conn = pymysql.connect(
            host=host,
            user=user,
            password=password,
            port=port,
            database=database,
            cursorclass=pymysql.cursors.DictCursor
        )

        with conn:
            with conn.cursor() as cursor:
                cursor.execute(query)
                results = cursor.fetchall()
                logger.info(f"Query returned {len(results)} rows.")
                return results

    except pymysql.MySQLError as e:
        logger.error(f"MySQL error: {e}")
        return []
