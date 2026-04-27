import pandas as pd
from app.utils.helpers import transform_metrics, compute_basic_metrics

def test_transform_metrics_realistic():
    
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
    
    
    assert set(result.columns) == {"genome_id", "protein_coding_genes", "transcript_count", "exon_count"}
    
    
    assert result.loc[result["genome_id"] == "GCA_000001405.28", "protein_coding_genes"].iloc[0] == 19973
    assert result.loc[result["genome_id"] == "GCA_000001405.28", "exon_count"].iloc[0] == 1374523

    
    assert result.loc[result["genome_id"] == "GCA_000001635.9", "transcript_count"].iloc[0] == 85341

def test_transform_metrics_invalid():
    
    df = pd.DataFrame({"genome_id": ["GCA_000001405.28"], "foo": ["A"], "bar": [10]})
    try:
        transform_metrics(df)
        assert False, "Expected ValueError for missing columns"
    except ValueError as e:
        assert "Missing required columns" in str(e)

def test_basic_metrics_realistic():
    
    df = pd.DataFrame({
        "genome_id": ["GCA_000242695.1", "GCA_000242695.1", "GCA_000242695.1"],
        "metric_name": ["ncrna_genes", "pseudogenes", "protein_coding_genes"],
        "metric_value": [1054, 381, 28144]
    })
    stats = compute_basic_metrics(df)
    
    
    assert stats == {
        "total_rows": 3, 
        "unique_genomes": 1, 
        "unique_metrics": 3
    }
