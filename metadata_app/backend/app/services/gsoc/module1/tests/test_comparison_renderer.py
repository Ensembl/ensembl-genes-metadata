"""
test_comparison_renderer.py

Unit tests for comparison_renderer.py.
All tests use synthetic fixture data — no database connection required.
"""

# pylint: disable=redefined-outer-name

from pathlib import Path

import pandas as pd
import pytest

from metadata_app.backend.app.services.gsoc.module1.comparison_renderer import (
    _build_comparison_rows,
    _busco_cell_style,
    _filter_rows,
    _fmt,
    _render_data_rows,
    _render_header_row,
    _status_badge_label,
    _status_badge_style,
    render_comparison_html,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def multi_row_anno_wide() -> pd.DataFrame:
    """anno_wide with two annotation rows for the same GCA: live + abandoned."""
    return pd.DataFrame(
        [
            {
                "gca": "GCA_963455335.1",
                "scientific_name": "Homo sapiens",
                "common_name": "human",
                "gb_status": "live",
                "annotation_method": "full_genebuild",
                "genebuilder": "genebuild_team",
                "annotated_version": 1.0,
                "release_date": "2025-01-16",
                "date_status_update": "2025-01-16",
                "protein_busco": "C:94.3%[S:91.2%,D:3.1%],F:2.1%,M:3.6%,n:255",
                "protein_busco_lineage": "primates_odb10",
                "assembly_busco": "C:98.1%[S:97.0%,D:1.1%],F:0.9%,M:1.0%,n:255",
                "assembly_busco_lineage": "primates_odb10",
                "coding_genes": 20442,
                "annotation_source": "ensembl",
                "bioproject_id": "PRJNA31257",
            },
            {
                "gca": "GCA_963455335.1",
                "scientific_name": "Homo sapiens",
                "common_name": "human",
                "gb_status": "abandoned",
                "annotation_method": "helixer",
                "genebuilder": "helixer_team",
                "annotated_version": None,
                "release_date": None,
                "date_status_update": "2026-03-04",
                "protein_busco": None,
                "protein_busco_lineage": None,
                "assembly_busco": None,
                "assembly_busco_lineage": None,
                "coding_genes": None,
                "annotation_source": None,
                "bioproject_id": None,
            },
        ]
    )


@pytest.fixture()
def single_row_anno_wide() -> pd.DataFrame:
    """anno_wide with one annotation row for a GCA."""
    return pd.DataFrame(
        [
            {
                "gca": "GCA_000000001.1",
                "scientific_name": "Sparse species",
                "common_name": None,
                "gb_status": "live",
                "annotation_method": "anno",
                "genebuilder": None,
                "annotated_version": 1.0,
                "release_date": None,
                "date_status_update": "2024-01-01",
                "protein_busco": None,
                "protein_busco_lineage": None,
                "assembly_busco": None,
                "assembly_busco_lineage": None,
                "coding_genes": None,
                "annotation_source": None,
                "bioproject_id": None,
            }
        ]
    )


@pytest.fixture()
def archive_row_anno_wide() -> pd.DataFrame:
    """anno_wide where the only row for a GCA has gb_status=archive."""
    return pd.DataFrame(
        [
            {
                "gca": "GCA_999999999.1",
                "scientific_name": "Old species",
                "common_name": None,
                "gb_status": "archive",
                "annotation_method": "full_genebuild",
                "genebuilder": None,
                "annotated_version": 1.0,
                "release_date": None,
                "date_status_update": "2020-01-01",
                "protein_busco": None,
                "protein_busco_lineage": None,
                "assembly_busco": None,
                "assembly_busco_lineage": None,
                "coding_genes": None,
                "annotation_source": None,
                "bioproject_id": None,
            }
        ]
    )


# ---------------------------------------------------------------------------
# _fmt
# ---------------------------------------------------------------------------


def test_fmt_none_returns_na() -> None:
    """None input returns N/A."""
    assert _fmt(None) == "N/A"


def test_fmt_empty_string_returns_na() -> None:
    """Empty string returns N/A."""
    assert _fmt("") == "N/A"


def test_fmt_none_string_returns_na() -> None:
    """The string 'none' returns N/A."""
    assert _fmt("none") == "N/A"


def test_fmt_nan_returns_na() -> None:
    """Float NaN returns N/A."""
    assert _fmt(float("nan")) == "N/A"


def test_fmt_normal_string() -> None:
    """A regular string is returned unchanged."""
    assert _fmt("full_genebuild") == "full_genebuild"


def test_fmt_integer() -> None:
    """An integer is cast to string."""
    assert _fmt(42) == "42"


# ---------------------------------------------------------------------------
# _status_badge_style and _status_badge_label
# ---------------------------------------------------------------------------


def test_badge_style_live() -> None:
    """Live status returns green badge style."""
    style = _status_badge_style("live")
    assert "166534" in style


def test_badge_style_abandoned() -> None:
    """Abandoned status returns grey badge style."""
    style = _status_badge_style("abandoned")
    assert "374151" in style


def test_badge_style_unknown_returns_default() -> None:
    """An unrecognised status returns the default grey style."""
    style = _status_badge_style("totally_unknown_status")
    assert "374151" in style


def test_badge_label_live() -> None:
    """Live status returns the Live label."""
    assert _status_badge_label("live") == "Live"


def test_badge_label_abandoned() -> None:
    """Abandoned status returns the Abandoned label."""
    assert _status_badge_label("abandoned") == "Abandoned"


def test_badge_label_unknown_returns_default() -> None:
    """An unrecognised status returns Unknown."""
    assert _status_badge_label("totally_unknown_status") == "Unknown"


# ---------------------------------------------------------------------------
# _busco_cell_style
# ---------------------------------------------------------------------------


def test_busco_cell_style_excellent() -> None:
    """A high BUSCO score returns a green background style."""
    style = _busco_cell_style(
        "C:97.0%[S:95.0%,D:2.0%],F:1.0%,M:2.0%,n:255", "protein_busco"
    )
    assert "dcfce7" in style


def test_busco_cell_style_poor() -> None:
    """A low BUSCO score returns an orange background style."""
    style = _busco_cell_style(
        "C:55.0%[S:50.0%,D:5.0%],F:5.0%,M:40.0%,n:255", "protein_busco"
    )
    assert "ffedd5" in style


def test_busco_cell_style_lineage_col_returns_empty() -> None:
    """Lineage columns are not colour-coded."""
    style = _busco_cell_style("primates_odb10", "protein_busco_lineage")
    assert style == ""


def test_busco_cell_style_non_busco_col_returns_empty() -> None:
    """Non-BUSCO columns return empty style."""
    style = _busco_cell_style("full_genebuild", "annotation_method")
    assert style == ""


def test_busco_cell_style_none_value_returns_empty() -> None:
    """None BUSCO value returns empty style."""
    style = _busco_cell_style(None, "protein_busco")
    assert style == ""


# ---------------------------------------------------------------------------
# _filter_rows
# ---------------------------------------------------------------------------


def test_filter_rows_returns_all_non_archive(multi_row_anno_wide: pd.DataFrame) -> None:
    """Both live and abandoned rows are returned; neither is archive."""
    result = _filter_rows(multi_row_anno_wide, "GCA_963455335.1")
    assert len(result) == 2


def test_filter_rows_excludes_archive(archive_row_anno_wide: pd.DataFrame) -> None:
    """Archive-only GCA returns an empty DataFrame."""
    result = _filter_rows(archive_row_anno_wide, "GCA_999999999.1")
    assert result.empty


def test_filter_rows_unknown_gca_returns_empty(
    multi_row_anno_wide: pd.DataFrame,
) -> None:
    """An unknown GCA returns an empty DataFrame."""
    result = _filter_rows(multi_row_anno_wide, "GCA_NOTEXIST.1")
    assert result.empty


def test_filter_rows_single_row(single_row_anno_wide: pd.DataFrame) -> None:
    """A GCA with one non-archive row returns a single-row DataFrame."""
    result = _filter_rows(single_row_anno_wide, "GCA_000000001.1")
    assert len(result) == 1


def test_filter_rows_resets_index(multi_row_anno_wide: pd.DataFrame) -> None:
    """The returned DataFrame has a clean 0-based index."""
    result = _filter_rows(multi_row_anno_wide, "GCA_963455335.1")
    assert list(result.index) == [0, 1]


# ---------------------------------------------------------------------------
# _build_comparison_rows
# ---------------------------------------------------------------------------


def test_build_comparison_rows_length(multi_row_anno_wide: pd.DataFrame) -> None:
    """Returns one tuple per column in _COMPARISON_COLUMNS."""
    rows = _filter_rows(multi_row_anno_wide, "GCA_963455335.1")
    result = _build_comparison_rows(rows)
    assert len(result) == 13


def test_build_comparison_rows_values_per_row(
    multi_row_anno_wide: pd.DataFrame,
) -> None:
    """Each tuple has one value per annotation row."""
    rows = _filter_rows(multi_row_anno_wide, "GCA_963455335.1")
    result = _build_comparison_rows(rows)
    for _, _, values in result:
        assert len(values) == 2


def test_build_comparison_rows_status_values(multi_row_anno_wide: pd.DataFrame) -> None:
    """The Status row contains live and abandoned values."""
    rows = _filter_rows(multi_row_anno_wide, "GCA_963455335.1")
    result = _build_comparison_rows(rows)
    status_row = next(r for r in result if r[0] == "Status")
    assert "live" in status_row[2]
    assert "abandoned" in status_row[2]


# ---------------------------------------------------------------------------
# _render_header_row
# ---------------------------------------------------------------------------


def test_render_header_row_contains_live_badge(
    multi_row_anno_wide: pd.DataFrame,
) -> None:
    """Header row contains the Live badge label."""
    rows = _filter_rows(multi_row_anno_wide, "GCA_963455335.1")
    header = _render_header_row(rows)
    assert "Live" in header


def test_render_header_row_contains_method(multi_row_anno_wide: pd.DataFrame) -> None:
    """Header row contains the annotation method."""
    rows = _filter_rows(multi_row_anno_wide, "GCA_963455335.1")
    header = _render_header_row(rows)
    assert "full_genebuild" in header
    assert "helixer" in header


def test_render_header_row_is_html(multi_row_anno_wide: pd.DataFrame) -> None:
    """Header row is wrapped in a tr tag."""
    rows = _filter_rows(multi_row_anno_wide, "GCA_963455335.1")
    header = _render_header_row(rows)
    assert header.startswith("<tr>")
    assert header.endswith("</tr>")


# ---------------------------------------------------------------------------
# _render_data_rows
# ---------------------------------------------------------------------------


def test_render_data_rows_contains_busco(multi_row_anno_wide: pd.DataFrame) -> None:
    """Data rows contain a BUSCO value from the live annotation."""
    rows = _filter_rows(multi_row_anno_wide, "GCA_963455335.1")
    comparison_rows = _build_comparison_rows(rows)
    html_out = _render_data_rows(comparison_rows)
    assert "94.3" in html_out


def test_render_data_rows_contains_na(multi_row_anno_wide: pd.DataFrame) -> None:
    """Data rows contain N/A for missing values in the abandoned row."""
    rows = _filter_rows(multi_row_anno_wide, "GCA_963455335.1")
    comparison_rows = _build_comparison_rows(rows)
    html_out = _render_data_rows(comparison_rows)
    assert "N/A" in html_out


# ---------------------------------------------------------------------------
# render_comparison_html — file creation
# ---------------------------------------------------------------------------


def test_render_comparison_html_creates_file(
    tmp_path: Path, multi_row_anno_wide: pd.DataFrame
) -> None:
    """render_comparison_html writes an HTML file to the output directory."""
    path = render_comparison_html("GCA_963455335.1", multi_row_anno_wide, tmp_path)
    assert path.exists()


def test_render_comparison_html_filename(
    tmp_path: Path, multi_row_anno_wide: pd.DataFrame
) -> None:
    """Output filename replaces dots with underscores and ends in _comparison.html."""
    path = render_comparison_html("GCA_963455335.1", multi_row_anno_wide, tmp_path)
    assert path.name == "GCA_963455335_1_comparison.html"


def test_render_comparison_html_returns_path(
    tmp_path: Path, multi_row_anno_wide: pd.DataFrame
) -> None:
    """render_comparison_html returns a Path object."""
    path = render_comparison_html("GCA_963455335.1", multi_row_anno_wide, tmp_path)
    assert isinstance(path, Path)


# ---------------------------------------------------------------------------
# render_comparison_html — HTML content
# ---------------------------------------------------------------------------


def test_render_comparison_html_contains_gca(
    tmp_path: Path, multi_row_anno_wide: pd.DataFrame
) -> None:
    """The HTML output contains the GCA accession."""
    render_comparison_html("GCA_963455335.1", multi_row_anno_wide, tmp_path)
    content = (tmp_path / "GCA_963455335_1_comparison.html").read_text(encoding="utf-8")
    assert "GCA_963455335.1" in content


def test_render_comparison_html_contains_species(
    tmp_path: Path, multi_row_anno_wide: pd.DataFrame
) -> None:
    """The HTML output contains the scientific name."""
    render_comparison_html("GCA_963455335.1", multi_row_anno_wide, tmp_path)
    content = (tmp_path / "GCA_963455335_1_comparison.html").read_text(encoding="utf-8")
    assert "Homo sapiens" in content


def test_render_comparison_html_contains_both_methods(
    tmp_path: Path, multi_row_anno_wide: pd.DataFrame
) -> None:
    """The HTML output contains both annotation methods."""
    render_comparison_html("GCA_963455335.1", multi_row_anno_wide, tmp_path)
    content = (tmp_path / "GCA_963455335_1_comparison.html").read_text(encoding="utf-8")
    assert "full_genebuild" in content
    assert "helixer" in content


def test_render_comparison_html_is_valid_skeleton(
    tmp_path: Path, multi_row_anno_wide: pd.DataFrame
) -> None:
    """The generated file is a complete HTML document."""
    render_comparison_html("GCA_963455335.1", multi_row_anno_wide, tmp_path)
    content = (tmp_path / "GCA_963455335_1_comparison.html").read_text(encoding="utf-8")
    assert content.startswith("<!DOCTYPE html>")
    assert "</html>" in content


def test_render_comparison_html_contains_footer(
    tmp_path: Path, multi_row_anno_wide: pd.DataFrame
) -> None:
    """The HTML output contains the Ensembl Genebuild Metadata footer."""
    render_comparison_html("GCA_963455335.1", multi_row_anno_wide, tmp_path)
    content = (tmp_path / "GCA_963455335_1_comparison.html").read_text(encoding="utf-8")
    assert "Ensembl Genebuild Metadata" in content


def test_render_comparison_html_single_row(
    tmp_path: Path, single_row_anno_wide: pd.DataFrame
) -> None:
    """render_comparison_html works for a GCA with only one annotation row."""
    path = render_comparison_html("GCA_000000001.1", single_row_anno_wide, tmp_path)
    assert path.exists()


# ---------------------------------------------------------------------------
# render_comparison_html — error handling
# ---------------------------------------------------------------------------


def test_render_comparison_html_raises_on_unknown_gca(
    tmp_path: Path, multi_row_anno_wide: pd.DataFrame
) -> None:
    """ValueError is raised when the GCA is not found in anno_wide."""
    with pytest.raises(ValueError, match="GCA_NOTEXIST"):
        render_comparison_html("GCA_NOTEXIST.1", multi_row_anno_wide, tmp_path)


def test_render_comparison_html_raises_on_archive_only(
    tmp_path: Path, archive_row_anno_wide: pd.DataFrame
) -> None:
    """ValueError is raised when all rows for the GCA are archived."""
    with pytest.raises(ValueError):
        render_comparison_html("GCA_999999999.1", archive_row_anno_wide, tmp_path)
