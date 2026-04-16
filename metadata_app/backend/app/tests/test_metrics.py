import pandas as pd
from app.utils.helpers import transform_metrics, compute_basic_metrics

def test_transform_metrics_basic():
    data = {
        "genome_id": [1, 1, 2, 2],
        "metric_name": ["A", "B", "A", "B"],
        "metric_value": [10, 20, 30, 40]
    }
    df = pd.DataFrame(data)
    result = transform_metrics(df)
    assert set(result.columns) == {"genome_id", "A", "B"}
    assert result.loc[result["genome_id"] == 1, "A"].iloc[0] == 10
    assert result.loc[result["genome_id"] == 2, "B"].iloc[0] == 40

def test_transform_metrics_invalid():
    df = pd.DataFrame({"genome_id": [1], "foo": ["A"], "bar": [10]})
    try:
        transform_metrics(df)
        assert False, "Expected ValueError for missing columns"
    except ValueError as e:
        assert "Missing required columns" in str(e)

def test_basic_metrics():
    df = pd.DataFrame({
        "genome_id": [1, 1, 2],
        "metric_name": ["A", "B", "A"],
        "metric_value": [10, 20, 30]
    })
    stats = compute_basic_metrics(df)
    assert stats == {"total_rows": 3, "unique_genomes": 2, "unique_metrics": 2}
