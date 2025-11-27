import pandas as pd
from logger_settings import get_logger
import pymysql

logger = get_logger(__name__)

def add_missing_methods(gb_status: pd.DataFrame, production_status: pd.DataFrame, test: bool, password: str) -> pd.DataFrame:
    """
    Fill missing annotation_method in live or handed_over rows using production data,
    and optionally update the database. Returns only the rows that were filled.
    """
    merge_keys = ["gca_accession", "release_date"]

    # Merge registry and production data
    df = pd.merge(
        gb_status,
        production_status,
        on=merge_keys,
        how='left',
        suffixes=('_registry', '_production')
    )

    # Fill missing annotation_method for live/handed_over rows
    condition_missing = (
        df['gb_status'].isin(['live', 'handed_over']) &
        df['annotation_method_registry'].isna() &
        df['annotation_method_production'].notna()
    )
    df.loc[condition_missing, 'annotation_method_new'] = df.loc[
        condition_missing, 'annotation_method_production'
    ]

    updated_df = df[condition_missing].copy()
    logger.info(f"Filled missing annotation_method for {len(updated_df)} rows")

    if not test and not updated_df.empty:
        try:
            connection = pymysql.connect(
                host="mysql-ens-genebuild-prod-1",
                user="ensadmin",
                password=password,
                port=4527,
                database="gb_assembly_metadata",
                autocommit=True
            )

            with connection.cursor() as cursor:
                for _, row in updated_df.iterrows():
                    sql = """
                        UPDATE genebuild_status
                        SET annotation_method = %s
                        WHERE genebuild_status_id = %s
                    """
                    cursor.execute(sql, (row['annotation_method_new'], row['genebuild_status_id']))

            logger.info(f"Updated {len(updated_df)} genebuild_status rows in the database")

        except pymysql.Error as e:
            logger.error(f"MySQL error: {e}")
            raise

        finally:
            connection.close()
    else:
        logger.info("No method updates applied.")
        if not updated_df.empty:
            logger.info("Rows that would be updated (method)\n%s", updated_df.to_string(index=False))

    return updated_df