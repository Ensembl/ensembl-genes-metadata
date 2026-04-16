import pandas as pd
from app.utils.helpers import transform_metrics, compute_basic_metrics

def test_transform_metrics_realistic():
    # Using realistic Ensembl genome IDs and annotation metrics
    data = {
        "genome_id": [
            "GCA_000001405.28", "GCA_000001405.28", "GCA_000001405.28",
            "GCA_000001635.9", "GCA_000001635.9", "GCA_000001635.9"
        ],
        "metric_name": [
            "protein_coding_genes", "transcript_count", "exon_count",
            "protein_coding_genes", "transcript_count", "exon_count"
        ],
        "metric_value": [
            19973, 237467, 1374523,
            22287, 85341, 624513
        ]
    }
    df = pd.DataFrame(data)
    result = transform_metrics(df)
    
    # Check if wide format correctly has the metric names as columns
    assert set(result.columns) == {"genome_id", "protein_coding_genes", "transcript_count", "exon_count"}
    
    # Check if values pivoted accurately for human (GCA_000001405.28)
    assert result.loc[result["genome_id"] == "GCA_000001405.28", "protein_coding_genes"].iloc[0] == 19973
    assert result.loc[result["genome_id"] == "GCA_000001405.28", "exon_count"].iloc[0] == 1374523

    # Check if values pivoted accurately for mouse (GCA_000001635.9)
    assert result.loc[result["genome_id"] == "GCA_000001635.9", "transcript_count"].iloc[0] == 85341

def test_transform_metrics_invalid():
    # If the user feeds a table missing required columns
    df = pd.DataFrame({"genome_id": ["GCA_000001405.28"], "foo": ["A"], "bar": [10]})
    try:
        transform_metrics(df)
        assert False, "Expected ValueError for missing columns"
    except ValueError as e:
        assert "Missing required columns" in str(e)

def test_basic_metrics_realistic():
    # Using realistic dataset fragment
    df = pd.DataFrame({
        "genome_id": ["GCA_000242695.1", "GCA_000242695.1", "GCA_000242695.1"],
        "metric_name": ["ncrna_genes", "pseudogenes", "protein_coding_genes"],
        "metric_value": [1054, 381, 28144]
    })
    stats = compute_basic_metrics(df)
    
    # Validating realistic counts 
    assert stats == {
        "total_rows": 3, 
        "unique_genomes": 1, 
        "unique_metrics": 3
    }
