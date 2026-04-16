import pandas as pd

def transform_metrics(df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert key-value metrics into structured tabular format.

    Expected columns:
    - genome_id
    - metric_name
    - metric_value
    """

    if df.empty:
        return df

    required_cols = {"genome_id", "metric_name", "metric_value"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"Missing required columns: {required_cols}")

    result = df.pivot(
        index="genome_id",
        columns="metric_name",
        values="metric_value"
    )

    return result.reset_index()

def compute_basic_metrics(df: pd.DataFrame) -> dict:
    """
    Compute total rows and unique count metrics.
    """
    if df.empty:
        return {}
    return {
        "total_rows": len(df),
        "unique_genomes": df["genome_id"].nunique(),
        "unique_metrics": df["metric_name"].nunique()
    }
