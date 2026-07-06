"""
test_genome_report.py

Unit tests for genome_report.py.
All tests use synthetic fixture data — no database connection required.
"""

# pylint: disable=redefined-outer-name

import pandas as pd
import pytest

from metadata_app.backend.app.services.gsoc.module1.genome_report import (
    GenomeReport,
    _busco_complete_as_float,
    _busco_extra_as_dict,
    _select_current_annotation_row,
    extract_genome_report,
    safe_float,
    safe_int,
    safe_str,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def full_row() -> dict:
    """A fully-populated anno_wide row for a mock genome."""
    return {
        "gca": "GCA_000001405.29",
        "scientific_name": "Homo sapiens",
        "common_name": "human",
        "lowest_taxon_id": 9606,
        "internal_clade": "Primates",
        "gb_status": "live",
        "annotation_method": "full_genebuild",
        "genebuilder": "genebuild_team",
        "annotation_source": "ensembl",
        "bioproject_id": "PRJNA31257",
        "associated_project": "Human Genome Project",
        "release_date": "2023-01-15",
        "date_status_update": "2023-01-15",
        "last_genebuild_update": "2023-01-10",
        "annotated_version": 29.0,
        "assembly_version": 29.0,
        "protein_busco": "C:94.3%[S:91.2%,D:3.1%],F:2.1%,M:3.6%,n:255",
        "protein_busco_lineage": "primates_odb10",
        "protein_busco_version": "5.4.3",
        "assembly_busco": "C:98.1%[S:97.0%,D:1.1%],F:0.9%,M:1.0%,n:255",
        "assembly_busco_lineage": "primates_odb10",
        "assembly_busco_version": "5.4.3",
        "coding_genes": 20442,
        "latest_annotated": "Yes",
        "ftp": "https://ftp.ebi.ac.uk/pub/ensemblorganisms/Homo_sapiens/GCA_000001405.29/",
    }


@pytest.fixture()
def sparse_row() -> dict:
    """A minimal anno_wide row with most optional fields missing."""
    return {
        "gca": "GCA_000000001.1",
        "scientific_name": "Sparse species",
        "common_name": None,
        "lowest_taxon_id": None,
        "internal_clade": None,
        "gb_status": None,
        "annotation_method": None,
        "genebuilder": None,
        "annotation_source": None,
        "bioproject_id": None,
        "associated_project": None,
        "release_date": None,
        "date_status_update": None,
        "last_genebuild_update": None,
        "annotated_version": None,
        "assembly_version": None,
        "protein_busco": None,
        "protein_busco_lineage": None,
        "protein_busco_version": None,
        "assembly_busco": None,
        "assembly_busco_lineage": None,
        "assembly_busco_version": None,
        "coding_genes": None,
        "latest_annotated": None,
        "ftp": None,
    }


@pytest.fixture()
def full_anno_wide(full_row: dict) -> pd.DataFrame:
    """A single-row anno_wide DataFrame for the full genome."""
    return pd.DataFrame([full_row])


@pytest.fixture()
def sparse_anno_wide(sparse_row: dict) -> pd.DataFrame:
    """A single-row anno_wide DataFrame for the sparse genome."""
    return pd.DataFrame([sparse_row])


# ---------------------------------------------------------------------------
# safe_str
# ---------------------------------------------------------------------------


def test_safe_str_none_returns_none() -> None:
    """None input returns None."""
    assert safe_str(None) is None


def test_safe_str_empty_string_returns_none() -> None:
    """Empty string input returns None."""
    assert safe_str("") is None


def test_safe_str_nan_returns_none() -> None:
    """NaN input returns None."""
    assert safe_str(float("nan")) is None


def test_safe_str_normal_string() -> None:
    """A non-empty string is returned unchanged."""
    assert safe_str("primates_odb10") == "primates_odb10"


def test_safe_str_int_cast() -> None:
    """An integer is cast to its string representation."""
    assert safe_str(42) == "42"


# ---------------------------------------------------------------------------
# safe_int
# ---------------------------------------------------------------------------


def test_safe_int_none_returns_none() -> None:
    """None input returns None."""
    assert safe_int(None) is None


def test_safe_int_empty_string_returns_none() -> None:
    """Empty string input returns None."""
    assert safe_int("") is None


def test_safe_int_valid_int() -> None:
    """A valid integer is returned as-is."""
    assert safe_int(9606) == 9606


def test_safe_int_float_truncates() -> None:
    """A float is truncated to integer."""
    assert safe_int(9606.9) == 9606


def test_safe_int_invalid_string_returns_none() -> None:
    """A non-numeric string returns None without raising."""
    assert safe_int("not_a_number") is None


# ---------------------------------------------------------------------------
# safe_float
# ---------------------------------------------------------------------------


def test_safe_float_none_returns_none() -> None:
    """None input returns None."""
    assert safe_float(None) is None


def test_safe_float_empty_string_returns_none() -> None:
    """Empty string input returns None."""
    assert safe_float("") is None


def test_safe_float_valid_float() -> None:
    """A valid float is returned as-is."""
    assert safe_float(94.3) == pytest.approx(94.3)


def test_safe_float_string_number() -> None:
    """A numeric string is correctly cast to float."""
    assert safe_float("94.3") == pytest.approx(94.3)


def test_safe_float_invalid_string_returns_none() -> None:
    """A non-numeric string returns None without raising."""
    assert safe_float("bad") is None


# ---------------------------------------------------------------------------
# _busco_complete_as_float
# ---------------------------------------------------------------------------


def test_busco_complete_as_float_returns_float() -> None:
    """A numeric complete value is returned as float."""
    parsed: dict = {"complete": 94.3, "extra": {}}
    assert _busco_complete_as_float(parsed) == pytest.approx(94.3)  # type: ignore[arg-type]


def test_busco_complete_as_float_none_returns_none() -> None:
    """None complete value returns None."""
    parsed: dict = {"complete": None, "extra": {}}
    assert _busco_complete_as_float(parsed) is None  # type: ignore[arg-type]


def test_busco_complete_as_float_dict_returns_none() -> None:
    """A dict value (the extra field type) must not be cast to float."""
    parsed: dict = {"complete": {"S": 91.2}, "extra": {}}
    assert _busco_complete_as_float(parsed) is None  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# _busco_extra_as_dict
# ---------------------------------------------------------------------------


def test_busco_extra_as_dict_returns_dict() -> None:
    """A dict extra value is returned as-is."""
    parsed: dict = {"complete": 94.3, "extra": {"E": 4.5}}
    assert _busco_extra_as_dict(parsed) == {"E": 4.5}  # type: ignore[arg-type]


def test_busco_extra_as_dict_none_returns_empty() -> None:
    """None extra value returns an empty dict."""
    parsed: dict = {"complete": 94.3, "extra": None}
    assert _busco_extra_as_dict(parsed) == {}  # type: ignore[arg-type]


def test_busco_extra_as_dict_non_dict_returns_empty() -> None:
    """A non-dict extra value returns an empty dict."""
    parsed: dict = {"complete": 94.3, "extra": 0.0}
    assert _busco_extra_as_dict(parsed) == {}  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# _select_current_annotation_row
# ---------------------------------------------------------------------------


def test_select_single_row_returns_it() -> None:
    """With one row, returns it directly without logging a warning."""
    df = pd.DataFrame([{"gca": "GCA_000001405.29", "date_status_update": "2023-01-15"}])
    result = _select_current_annotation_row(df, "GCA_000001405.29")
    assert result["gca"] == "GCA_000001405.29"


def test_select_most_recent_row() -> None:
    """With two rows, the one with the later date_status_update is selected."""
    df = pd.DataFrame(
        [
            {
                "gca": "GCA_000001405.29",
                "date_status_update": "2022-01-01",
                "gb_status": "completed",
            },
            {
                "gca": "GCA_000001405.29",
                "date_status_update": "2023-06-15",
                "gb_status": "live",
            },
        ]
    )
    result = _select_current_annotation_row(df, "GCA_000001405.29")
    assert result["gb_status"] == "live"


def test_select_falls_back_on_missing_date_column() -> None:
    """Falls back to first row when date_status_update column is absent."""
    df = pd.DataFrame(
        [
            {"gca": "GCA_000001405.29", "gb_status": "completed"},
            {"gca": "GCA_000001405.29", "gb_status": "live"},
        ]
    )
    result = _select_current_annotation_row(df, "GCA_000001405.29")
    assert result["gb_status"] == "completed"


def test_select_falls_back_on_all_unparseable_dates() -> None:
    """Falls back to first row when no date_status_update value is parseable."""
    df = pd.DataFrame(
        [
            {
                "gca": "GCA_000001405.29",
                "date_status_update": "not-a-date",
                "gb_status": "completed",
            },
            {
                "gca": "GCA_000001405.29",
                "date_status_update": "also-bad",
                "gb_status": "live",
            },
        ]
    )
    result = _select_current_annotation_row(df, "GCA_000001405.29")
    assert result["gb_status"] == "completed"


# ---------------------------------------------------------------------------
# extract_genome_report
# ---------------------------------------------------------------------------


def test_extract_genome_report_returns_dataclass(full_anno_wide: pd.DataFrame) -> None:
    """extract_genome_report returns a GenomeReport instance."""
    report = extract_genome_report("GCA_000001405.29", full_anno_wide)
    assert isinstance(report, GenomeReport)


def test_extract_genome_report_gca(full_anno_wide: pd.DataFrame) -> None:
    """The GCA accession is copied to the report correctly."""
    report = extract_genome_report("GCA_000001405.29", full_anno_wide)
    assert report.gca == "GCA_000001405.29"


def test_extract_genome_report_scientific_name(full_anno_wide: pd.DataFrame) -> None:
    """The scientific name is extracted correctly from the anno_wide row."""
    report = extract_genome_report("GCA_000001405.29", full_anno_wide)
    assert report.scientific_name == "Homo sapiens"


def test_extract_genome_report_protein_busco_complete(
    full_anno_wide: pd.DataFrame,
) -> None:
    """Protein BUSCO completeness is extracted and parsed correctly."""
    report = extract_genome_report("GCA_000001405.29", full_anno_wide)
    assert report.protein_busco_complete == pytest.approx(94.3)


def test_extract_genome_report_assembly_busco_complete(
    full_anno_wide: pd.DataFrame,
) -> None:
    """Assembly BUSCO completeness is extracted correctly from the anno_wide row."""
    report = extract_genome_report("GCA_000001405.29", full_anno_wide)
    assert report.assembly_busco_complete == pytest.approx(98.1)


def test_extract_genome_report_coding_genes(full_anno_wide: pd.DataFrame) -> None:
    """Coding gene count is extracted correctly from the anno_wide row."""
    report = extract_genome_report("GCA_000001405.29", full_anno_wide)
    assert report.coding_genes == 20442


def test_extract_genome_report_quality_label(full_anno_wide: pd.DataFrame) -> None:
    """A 94.3% BUSCO score should map to the Good quality band."""
    report = extract_genome_report("GCA_000001405.29", full_anno_wide)
    assert report.protein_busco_quality == "Good"


def test_extract_genome_report_raises_on_missing_gca(
    full_anno_wide: pd.DataFrame,
) -> None:
    """extract_genome_report raises ValueError for an unknown GCA."""
    with pytest.raises(ValueError, match="GCA_MISSING"):
        extract_genome_report("GCA_MISSING", full_anno_wide)


def test_extract_genome_report_sparse_no_crash(sparse_anno_wide: pd.DataFrame) -> None:
    """extract_genome_report handles a row with all-None optional fields."""
    report = extract_genome_report("GCA_000000001.1", sparse_anno_wide)
    assert report.gca == "GCA_000000001.1"
    assert report.protein_busco_complete is None
    assert report.coding_genes is None


def test_extract_genome_report_busco_extra_populated(
    full_anno_wide: pd.DataFrame,
) -> None:
    """protein_busco_extra is a dict (empty when no extra fields in string)."""
    report = extract_genome_report("GCA_000001405.29", full_anno_wide)
    assert isinstance(report.protein_busco_extra, dict)


def test_extract_genome_report_unknown_quality_when_no_busco(
    sparse_anno_wide: pd.DataFrame,
) -> None:
    """When protein_busco is None, quality label is Unknown."""
    report = extract_genome_report("GCA_000000001.1", sparse_anno_wide)
    assert report.protein_busco_quality == "Unknown"
