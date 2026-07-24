"""
test_clade_analysis.py

Unit tests for clade_analysis.py.
All tests use synthetic fixture data -- no database connection required.
"""

# pylint: disable=redefined-outer-name

import numpy as np
import pandas as pd
import pytest

from metadata_app.backend.app.services.gsoc.module2.clade_analysis import (
    OutlierResult,
    _busco_to_pct,
    _mad_scores,
    _prepare_features,
    _run_pca,
    outlier_results_to_dataframe,
    run_clade_analysis,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def minimal_clade_df() -> pd.DataFrame:
    """A small synthetic clade metrics DataFrame with two clades."""
    rows = []
    for i in range(15):
        rows.append(
            {
                "gca": f"GCA_mammal_{i}.1",
                "scientific_name": f"Mammal species {i}",
                "lowest_taxon_id": 9000 + i,
                "annotation_method": "full_genebuild",
                "clade": "mammalia",
                "genebuild.busco_completeness": str(13000 + i * 10),
                "genebuild.busco_duplicated": str(200 + i),
                "genebuild.busco_fragmented": str(100 + i),
                "genebuild.busco_missing": str(400 + i),
                "genebuild.busco_single_copy": str(12800 + i * 10),
                "genebuild.stats.coding_genes": str(20000 + i * 100),
                "genebuild.stats.total_transcripts": str(50000 + i * 200),
                "genebuild.stats.transcripts_per_gene": "2.5",
                "genebuild.stats.average_cds_length": "1500",
                "genebuild.stats.average_coding_intron_length": "5000",
                "genebuild.stats.single_exon_coding_genes": "2000",
                "genebuild.stats.overlapping_coding_genes": "500",
            }
        )
    for i in range(12):
        rows.append(
            {
                "gca": f"GCA_aves_{i}.1",
                "scientific_name": f"Bird species {i}",
                "lowest_taxon_id": 8000 + i,
                "annotation_method": "full_genebuild",
                "clade": "aves",
                "genebuild.busco_completeness": str(11000 + i * 10),
                "genebuild.busco_duplicated": str(150 + i),
                "genebuild.busco_fragmented": str(80 + i),
                "genebuild.busco_missing": str(300 + i),
                "genebuild.busco_single_copy": str(10850 + i * 10),
                "genebuild.stats.coding_genes": str(15000 + i * 100),
                "genebuild.stats.total_transcripts": str(30000 + i * 200),
                "genebuild.stats.transcripts_per_gene": "2.0",
                "genebuild.stats.average_cds_length": "1200",
                "genebuild.stats.average_coding_intron_length": "3000",
                "genebuild.stats.single_exon_coding_genes": "1500",
                "genebuild.stats.overlapping_coding_genes": "300",
            }
        )
    return pd.DataFrame(rows)


@pytest.fixture()
def outlier_clade_df() -> pd.DataFrame:
    """A DataFrame where one genome is a clear outlier in mammalia."""
    rows = []
    for i in range(14):
        rows.append(
            {
                "gca": f"GCA_normal_{i}.1",
                "scientific_name": f"Normal mammal {i}",
                "lowest_taxon_id": 9000 + i,
                "annotation_method": "full_genebuild",
                "clade": "mammalia",
                "genebuild.busco_completeness": str(13000 + i * 5),
                "genebuild.busco_duplicated": str(200 + i * 2),
                "genebuild.busco_fragmented": str(100 + i),
                "genebuild.busco_missing": str(400 + i * 3),
                "genebuild.busco_single_copy": str(12800 + i * 5),
                "genebuild.stats.coding_genes": str(20000 + i * 50),
                "genebuild.stats.total_transcripts": str(50000 + i * 100),
                "genebuild.stats.transcripts_per_gene": str(2.5 + i * 0.01),
                "genebuild.stats.average_cds_length": str(1500 + i * 10),
                "genebuild.stats.average_coding_intron_length": str(5000 + i * 20),
                "genebuild.stats.single_exon_coding_genes": str(2000 + i * 5),
                "genebuild.stats.overlapping_coding_genes": str(500 + i * 2),
            }
        )
    # One extreme outlier with very low BUSCO
    rows.append(
        {
            "gca": "GCA_outlier.1",
            "scientific_name": "Outlier mammal",
            "lowest_taxon_id": 9999,
            "annotation_method": "full_genebuild",
            "clade": "mammalia",
            "genebuild.busco_completeness": "1000",
            "genebuild.busco_duplicated": "50",
            "genebuild.busco_fragmented": "500",
            "genebuild.busco_missing": "5000",
            "genebuild.busco_single_copy": "950",
            "genebuild.stats.coding_genes": "5000",
            "genebuild.stats.total_transcripts": "8000",
            "genebuild.stats.transcripts_per_gene": "1.2",
            "genebuild.stats.average_cds_length": "800",
            "genebuild.stats.average_coding_intron_length": "1000",
            "genebuild.stats.single_exon_coding_genes": "3000",
            "genebuild.stats.overlapping_coding_genes": "100",
        }
    )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# _busco_to_pct
# ---------------------------------------------------------------------------


def test_busco_to_pct_adds_columns(minimal_clade_df: pd.DataFrame) -> None:
    """_busco_to_pct adds the four percentage columns."""
    result = _busco_to_pct(minimal_clade_df)
    for col in [
        "busco_completeness_pct",
        "busco_duplicated_pct",
        "busco_fragmented_pct",
        "busco_missing_pct",
    ]:
        assert col in result.columns


def test_busco_to_pct_range(minimal_clade_df: pd.DataFrame) -> None:
    """All percentage values should be between 0 and 100."""
    result = _busco_to_pct(minimal_clade_df)
    assert (result["busco_completeness_pct"].dropna() >= 0).all()
    assert (result["busco_completeness_pct"].dropna() <= 100).all()


def test_busco_to_pct_sums_to_100(minimal_clade_df: pd.DataFrame) -> None:
    """Completeness + fragmented + missing should sum to ~100."""
    result = _busco_to_pct(minimal_clade_df)
    total = (
        result["busco_completeness_pct"]
        + result["busco_fragmented_pct"]
        + result["busco_missing_pct"]
    )
    assert (total.dropna() - 100).abs().max() < 0.01


def test_busco_to_pct_does_not_modify_input(minimal_clade_df: pd.DataFrame) -> None:
    """_busco_to_pct does not modify the input DataFrame."""
    original_cols = list(minimal_clade_df.columns)
    _busco_to_pct(minimal_clade_df)
    assert list(minimal_clade_df.columns) == original_cols


# ---------------------------------------------------------------------------
# _mad_scores
# ---------------------------------------------------------------------------


def test_mad_scores_returns_array() -> None:
    """_mad_scores returns a numpy array."""
    values = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    result = _mad_scores(values)
    assert isinstance(result, np.ndarray)


def test_mad_scores_median_is_zero() -> None:
    """The median value always has a MAD score of 0."""
    values = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    result = _mad_scores(values)
    assert result[2] == pytest.approx(0.0)


def test_mad_scores_outlier_has_high_score() -> None:
    """An extreme value should have a much higher MAD score than others."""
    values = np.array([10.0, 10.1, 10.0, 9.9, 10.0, 10.1, 10.0, 9.9, 10.0, 100.0])
    result = _mad_scores(values)
    assert result[-1] > result[0]


def test_mad_scores_zero_mad_returns_zeros() -> None:
    """When all values are identical, MAD is 0 and scores are all 0."""
    values = np.array([5.0, 5.0, 5.0, 5.0, 5.0])
    result = _mad_scores(values)
    assert (result == 0.0).all()


def test_mad_scores_all_positive() -> None:
    """MAD scores are always non-negative."""
    values = np.array([1.0, 2.0, 3.0, 100.0, 2.5])
    result = _mad_scores(values)
    assert (result >= 0).all()


# ---------------------------------------------------------------------------
# _prepare_features
# ---------------------------------------------------------------------------


def test_prepare_features_converts_to_numeric(minimal_clade_df: pd.DataFrame) -> None:
    """_prepare_features converts string columns to numeric."""
    df = _busco_to_pct(minimal_clade_df)
    result = _prepare_features(df)
    assert result["genebuild.stats.coding_genes"].dtype in [np.float64, np.int64, float]


def test_prepare_features_fills_nulls(minimal_clade_df: pd.DataFrame) -> None:
    """_prepare_features fills null values with column medians."""
    df = _busco_to_pct(minimal_clade_df)
    df.loc[0, "genebuild.stats.coding_genes"] = None
    result = _prepare_features(df)
    assert result["genebuild.stats.coding_genes"].isna().sum() == 0


# ---------------------------------------------------------------------------
# _run_pca
# ---------------------------------------------------------------------------


def test_run_pca_returns_correct_shape() -> None:
    """PCA output has shape (n_samples, n_components)."""
    matrix = np.random.rand(20, 5)
    result = _run_pca(matrix, n_components=2)
    assert result.shape == (20, 2)


def test_run_pca_single_component() -> None:
    """PCA works with n_components=1."""
    matrix = np.random.rand(15, 4)
    result = _run_pca(matrix, n_components=1)
    assert result.shape == (15, 1)


# ---------------------------------------------------------------------------
# run_clade_analysis
# ---------------------------------------------------------------------------


def test_run_clade_analysis_returns_dict(minimal_clade_df: pd.DataFrame) -> None:
    """run_clade_analysis returns a dict."""
    result = run_clade_analysis(minimal_clade_df)
    assert isinstance(result, dict)


def test_run_clade_analysis_has_both_clades(minimal_clade_df: pd.DataFrame) -> None:
    """Both clades in the fixture appear in the results."""
    result = run_clade_analysis(minimal_clade_df)
    assert "mammalia" in result
    assert "aves" in result


def test_run_clade_analysis_returns_outlier_results(
    minimal_clade_df: pd.DataFrame,
) -> None:
    """Each clade maps to a list of OutlierResult objects."""
    result = run_clade_analysis(minimal_clade_df)
    for clade_results in result.values():
        assert isinstance(clade_results, list)
        assert all(isinstance(r, OutlierResult) for r in clade_results)


def test_run_clade_analysis_empty_input() -> None:
    """Empty DataFrame returns empty dict."""
    result = run_clade_analysis(pd.DataFrame())
    assert not result


def test_run_clade_analysis_flags_outlier(outlier_clade_df: pd.DataFrame) -> None:
    """The extreme outlier genome is flagged as an outlier."""
    result = run_clade_analysis(outlier_clade_df)
    mammalia_results = result.get("mammalia", [])
    outlier = next((r for r in mammalia_results if r.gca == "GCA_outlier.1"), None)
    assert outlier is not None
    assert outlier.is_outlier is True


def test_run_clade_analysis_normal_not_flagged(outlier_clade_df: pd.DataFrame) -> None:
    """Normal genomes are not flagged as outliers."""
    result = run_clade_analysis(outlier_clade_df)
    mammalia_results = result.get("mammalia", [])
    normal = [r for r in mammalia_results if r.gca != "GCA_outlier.1"]
    assert all(not r.is_outlier for r in normal)


def test_run_clade_analysis_result_count(minimal_clade_df: pd.DataFrame) -> None:
    """Number of results per clade matches number of genomes in that clade."""
    result = run_clade_analysis(minimal_clade_df)
    assert len(result["mammalia"]) == 15
    assert len(result["aves"]) == 12


# ---------------------------------------------------------------------------
# outlier_results_to_dataframe
# ---------------------------------------------------------------------------


def test_outlier_results_to_dataframe_returns_df(
    minimal_clade_df: pd.DataFrame,
) -> None:
    """outlier_results_to_dataframe returns a DataFrame."""
    results = run_clade_analysis(minimal_clade_df)
    df = outlier_results_to_dataframe(results)
    assert isinstance(df, pd.DataFrame)


def test_outlier_results_to_dataframe_has_expected_columns(
    minimal_clade_df: pd.DataFrame,
) -> None:
    """The output DataFrame has all expected columns."""
    results = run_clade_analysis(minimal_clade_df)
    df = outlier_results_to_dataframe(results)
    for col in ["gca", "scientific_name", "clade", "is_outlier", "mad_score"]:
        assert col in df.columns


def test_outlier_results_to_dataframe_empty_input() -> None:
    """Empty results dict returns empty DataFrame."""
    df = outlier_results_to_dataframe({})
    assert df.empty


def test_outlier_results_to_dataframe_row_count(
    minimal_clade_df: pd.DataFrame,
) -> None:
    """Total rows equals total genomes across all clades."""
    results = run_clade_analysis(minimal_clade_df)
    df = outlier_results_to_dataframe(results)
    assert len(df) == len(minimal_clade_df)
