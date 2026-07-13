"""
test_clade_loader.py

Unit tests for clade_loader.py.
All tests use synthetic fixture data -- no database connection required.
"""

# pylint: disable=redefined-outer-name

import pandas as pd
import pytest

from metadata_app.backend.app.services.gsoc.module2.clade_loader import (
    MIN_CLADE_SIZE,
    _assign_clades,
    _build_placeholders,
    _drop_small_clades,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def sample_clade_data() -> dict:
    """Minimal clade_settings.json structure for testing."""
    return {
        "mammalia": {"taxon_id": 40674},
        "aves": {"taxon_id": 8782},
        "primates": {"taxon_id": 1758},
    }


@pytest.fixture()
def sample_taxonomy_dict() -> dict:
    """Taxonomy hierarchy for a few test taxa."""
    return {
        "9606": [
            {"taxon_class": "species", "taxon_class_id": 9606},
            {"taxon_class": "genus", "taxon_class_id": 9605},
            {"taxon_class": "family", "taxon_class_id": 9604},
            {"taxon_class": "order", "taxon_class_id": 9443},
            {"taxon_class": "class", "taxon_class_id": 40674},
        ],
        "9598": [
            {"taxon_class": "species", "taxon_class_id": 9598},
            {"taxon_class": "genus", "taxon_class_id": 9596},
            {"taxon_class": "class", "taxon_class_id": 40674},
        ],
        "9031": [
            {"taxon_class": "species", "taxon_class_id": 9031},
            {"taxon_class": "class", "taxon_class_id": 8782},
        ],
    }


@pytest.fixture()
def sample_base_df() -> pd.DataFrame:
    """A small anno_wide style DataFrame with three genomes."""
    return pd.DataFrame(
        [
            {
                "gca": "GCA_000001405.29",
                "lowest_taxon_id": 9606,
                "scientific_name": "Homo sapiens",
                "annotation_method": "full_genebuild",
            },
            {
                "gca": "GCA_000001515.5",
                "lowest_taxon_id": 9598,
                "scientific_name": "Pan troglodytes",
                "annotation_method": "full_genebuild",
            },
            {
                "gca": "GCA_000002315.5",
                "lowest_taxon_id": 9031,
                "scientific_name": "Gallus gallus",
                "annotation_method": "full_genebuild",
            },
        ]
    )


# ---------------------------------------------------------------------------
# _build_placeholders
# ---------------------------------------------------------------------------


def test_build_placeholders_single() -> None:
    """Single metric produces a quoted string."""
    result = _build_placeholders(["genebuild.busco"])
    assert result == "'genebuild.busco'"


def test_build_placeholders_multiple() -> None:
    """Multiple metrics are comma-separated and quoted."""
    result = _build_placeholders(["genebuild.busco", "genebuild.busco_dataset"])
    assert result == "'genebuild.busco', 'genebuild.busco_dataset'"


def test_build_placeholders_empty() -> None:
    """Empty list returns empty string."""
    result = _build_placeholders([])
    assert result == ""


# ---------------------------------------------------------------------------
# _assign_clades
# ---------------------------------------------------------------------------


def test_assign_clades_mammalia(
    sample_base_df: pd.DataFrame,
    sample_taxonomy_dict: dict,
    sample_clade_data: dict,
) -> None:
    """Homo sapiens and Pan troglodytes should be assigned to mammalia."""
    result = _assign_clades(sample_base_df, sample_taxonomy_dict, sample_clade_data)
    human = result[result["gca"] == "GCA_000001405.29"].iloc[0]
    chimp = result[result["gca"] == "GCA_000001515.5"].iloc[0]
    assert human["clade"] == "mammalia"
    assert chimp["clade"] == "mammalia"


def test_assign_clades_aves(
    sample_base_df: pd.DataFrame,
    sample_taxonomy_dict: dict,
    sample_clade_data: dict,
) -> None:
    """Gallus gallus should be assigned to aves."""
    result = _assign_clades(sample_base_df, sample_taxonomy_dict, sample_clade_data)
    chicken = result[result["gca"] == "GCA_000002315.5"].iloc[0]
    assert chicken["clade"] == "aves"


def test_assign_clades_unassigned_when_no_match() -> None:
    """A genome with no matching clade gets Unassigned."""
    df = pd.DataFrame(
        [
            {
                "gca": "GCA_999999999.1",
                "lowest_taxon_id": 99999,
                "scientific_name": "Unknown species",
                "annotation_method": "anno",
            }
        ]
    )
    taxonomy_dict = {"99999": [{"taxon_class": "species", "taxon_class_id": 99999}]}
    clade_data = {"mammalia": {"taxon_id": 40674}}
    result = _assign_clades(df, taxonomy_dict, clade_data)
    assert result.iloc[0]["clade"] == "Unassigned"


def test_assign_clades_adds_clade_column(
    sample_base_df: pd.DataFrame,
    sample_taxonomy_dict: dict,
    sample_clade_data: dict,
) -> None:
    """_assign_clades adds a clade column to the DataFrame."""
    result = _assign_clades(sample_base_df, sample_taxonomy_dict, sample_clade_data)
    assert "clade" in result.columns


def test_assign_clades_does_not_modify_input(
    sample_base_df: pd.DataFrame,
    sample_taxonomy_dict: dict,
    sample_clade_data: dict,
) -> None:
    """_assign_clades does not modify the input DataFrame."""
    original_cols = list(sample_base_df.columns)
    _assign_clades(sample_base_df, sample_taxonomy_dict, sample_clade_data)
    assert list(sample_base_df.columns) == original_cols


# ---------------------------------------------------------------------------
# _drop_small_clades
# ---------------------------------------------------------------------------


def test_drop_small_clades_removes_unassigned() -> None:
    """Unassigned genomes are always removed."""
    df = pd.DataFrame(
        [{"gca": f"GCA_{i}.1", "clade": "mammalia"} for i in range(MIN_CLADE_SIZE)]
        + [{"gca": "GCA_999.1", "clade": "Unassigned"}]
    )
    result = _drop_small_clades(df)
    assert "Unassigned" not in result["clade"].values


def test_drop_small_clades_removes_small_clades() -> None:
    """Clades with fewer than MIN_CLADE_SIZE genomes are removed."""
    df = pd.DataFrame(
        [{"gca": f"GCA_{i}.1", "clade": "mammalia"} for i in range(MIN_CLADE_SIZE)]
        + [{"gca": "GCA_small.1", "clade": "tiny_clade"}]
    )
    result = _drop_small_clades(df)
    assert "tiny_clade" not in result["clade"].values


def test_drop_small_clades_keeps_large_clades() -> None:
    """Clades with at least MIN_CLADE_SIZE genomes are kept."""
    df = pd.DataFrame(
        [{"gca": f"GCA_{i}.1", "clade": "mammalia"} for i in range(MIN_CLADE_SIZE)]
    )
    result = _drop_small_clades(df)
    assert "mammalia" in result["clade"].values


def test_drop_small_clades_resets_index() -> None:
    """The returned DataFrame has a clean 0-based index."""
    df = pd.DataFrame(
        [{"gca": f"GCA_{i}.1", "clade": "mammalia"} for i in range(MIN_CLADE_SIZE)]
    )
    result = _drop_small_clades(df)
    assert list(result.index) == list(range(len(result)))


def test_drop_small_clades_empty_input() -> None:
    """Empty DataFrame input returns empty DataFrame."""
    df = pd.DataFrame(columns=["gca", "clade"])
    result = _drop_small_clades(df)
    assert result.empty
