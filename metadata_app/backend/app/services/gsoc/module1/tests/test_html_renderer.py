"""
test_html_renderer.py

Unit tests for html_renderer.py.
All tests use synthetic fixture data — no database connection required.
"""

# pylint: disable=redefined-outer-name

from pathlib import Path

import pytest

from metadata_app.backend.app.services.gsoc.module1.genome_report import GenomeReport
from metadata_app.backend.app.services.gsoc.module1.html_renderer import (
    _busco_bar_color,
    _busco_chart_js,
    _fmt,
    _quality_badge_style,
    _safe_float,
    render_html,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def full_report() -> GenomeReport:
    """A fully-populated GenomeReport for a mock primate genome."""
    return GenomeReport(
        gca="GCA_000001515.5",
        scientific_name="Pan troglodytes",
        common_name="chimpanzee",
        lowest_taxon_id=9598,
        internal_clade="Primates",
        annotation_method="full_genebuild",
        genebuilder="genebuild_team",
        gb_status="live",
        annotation_source="ensembl",
        bioproject_id="PRJNA10627",
        associated_project="Chimp Genome Project",
        release_date="2023-06-01",
        date_status_update="2023-06-01",
        last_genebuild_update="2023-05-28",
        protein_busco_raw="C:97.7%[S:95.1%,D:2.6%],F:0.9%,M:1.4%,n:255",
        protein_busco_lineage="primates_odb10",
        protein_busco_version="5.4.3",
        protein_busco_complete=97.7,
        protein_busco_quality="Excellent",
        protein_busco_extra={},
        assembly_busco_raw="C:98.5%[S:97.2%,D:1.3%],F:0.7%,M:0.8%,n:255",
        assembly_busco_lineage="primates_odb10",
        assembly_busco_version="5.4.3",
        assembly_busco_complete=98.5,
        assembly_busco_extra={},
        coding_genes=23534,
        latest_annotated="Yes",
        annotated_version=5.0,
        assembly_version=5.0,
        ftp="https://ftp.ebi.ac.uk/pub/ensemblorganisms/Pan_troglodytes/GCA_000001515.5/",
    )


@pytest.fixture()
def sparse_report() -> GenomeReport:
    """A minimal GenomeReport with most optional fields as None."""
    return GenomeReport(
        gca="GCA_000000001.1",
        scientific_name="Sparse species",
        common_name=None,
        lowest_taxon_id=None,
        internal_clade=None,
        annotation_method=None,
        genebuilder=None,
        gb_status=None,
        annotation_source=None,
        bioproject_id=None,
        associated_project=None,
        release_date=None,
        date_status_update=None,
        last_genebuild_update=None,
        protein_busco_raw=None,
        protein_busco_lineage=None,
        protein_busco_version=None,
        protein_busco_complete=None,
        protein_busco_quality="Very Low",
        protein_busco_extra={},
        assembly_busco_raw=None,
        assembly_busco_lineage=None,
        assembly_busco_version=None,
        assembly_busco_complete=None,
        assembly_busco_extra={},
        coding_genes=None,
        latest_annotated=None,
        annotated_version=None,
        assembly_version=None,
        ftp=None,
    )


# ---------------------------------------------------------------------------
# _fmt
# ---------------------------------------------------------------------------


def test_fmt_none_returns_na() -> None:
    """None input returns the string N/A."""
    assert _fmt(None) == "N/A"


def test_fmt_empty_string_returns_na() -> None:
    """Empty string input returns the string N/A."""
    assert _fmt("") == "N/A"


def test_fmt_int_returns_string() -> None:
    """Integer input is cast to its string representation."""
    assert _fmt(42) == "42"


def test_fmt_float_returns_string() -> None:
    """Float input is cast to its string representation."""
    assert _fmt(97.7) == "97.7"


def test_fmt_normal_string_unchanged() -> None:
    """A non-empty string is returned unchanged."""
    assert _fmt("primates_odb10") == "primates_odb10"


# ---------------------------------------------------------------------------
# _safe_float
# ---------------------------------------------------------------------------


def test_safe_float_none_returns_zero() -> None:
    """None input returns 0.0."""
    assert _safe_float(None) == 0.0


def test_safe_float_valid_int() -> None:
    """An integer input is returned as a float."""
    assert _safe_float(95) == 95.0


def test_safe_float_valid_string_number() -> None:
    """A numeric string is correctly cast to float."""
    assert _safe_float("3.14") == 3.14


def test_safe_float_invalid_string_returns_zero() -> None:
    """A non-numeric string returns 0.0 without raising."""
    assert _safe_float("not_a_number") == 0.0


# ---------------------------------------------------------------------------
# _busco_bar_color
# ---------------------------------------------------------------------------


def test_busco_bar_color_excellent() -> None:
    """A BUSCO score of 97 returns the green colour."""
    assert _busco_bar_color(97.0) == "#16a34a"


def test_busco_bar_color_high() -> None:
    """A BUSCO score of 88 returns the amber colour."""
    assert _busco_bar_color(88.0) == "#d97706"


def test_busco_bar_color_medium() -> None:
    """A BUSCO score of 75 returns the orange colour."""
    assert _busco_bar_color(75.0) == "#ea580c"


def test_busco_bar_color_low() -> None:
    """A BUSCO score of 60 returns the red colour."""
    assert _busco_bar_color(60.0) == "#dc2626"


def test_busco_bar_color_boundary_95() -> None:
    """Exactly 95.0 maps to the green (Excellent) colour."""
    assert _busco_bar_color(95.0) == "#16a34a"


def test_busco_bar_color_boundary_85() -> None:
    """Exactly 85.0 maps to the amber (High) colour."""
    assert _busco_bar_color(85.0) == "#d97706"


def test_busco_bar_color_boundary_70() -> None:
    """Exactly 70.0 maps to the orange (Medium) colour."""
    assert _busco_bar_color(70.0) == "#ea580c"


# ---------------------------------------------------------------------------
# _quality_badge_style
# ---------------------------------------------------------------------------


def test_badge_style_excellent_contains_green() -> None:
    """Excellent quality badge contains green colour values."""
    style = _quality_badge_style("Excellent")
    assert "#166534" in style
    assert "#dcfce7" in style


def test_badge_style_very_low_contains_red() -> None:
    """Very Low quality badge contains a red colour value."""
    style = _quality_badge_style("Very Low")
    assert "#7f1d1d" in style


def test_badge_style_unknown_quality_returns_fallback() -> None:
    """An unrecognised quality label returns the neutral fallback colours."""
    style = _quality_badge_style("Unknown")
    assert "#374151" in style
    assert "#f3f4f6" in style


def test_badge_style_format_is_css() -> None:
    """The returned string is a valid inline CSS snippet."""
    style = _quality_badge_style("High")
    assert "color:" in style
    assert "background:" in style


# ---------------------------------------------------------------------------
# _busco_chart_js
# ---------------------------------------------------------------------------


def test_busco_chart_js_contains_chart_call(full_report: GenomeReport) -> None:
    """The generated JS contains a Chart.js instantiation call."""
    js = _busco_chart_js(full_report)
    assert "new Chart" in js


def test_busco_chart_js_contains_single_value(full_report: GenomeReport) -> None:
    """The single-copy BUSCO percentage is present in the JS output."""
    js = _busco_chart_js(full_report)
    assert "95.1" in js


def test_busco_chart_js_four_datasets(full_report: GenomeReport) -> None:
    """The chart script contains exactly four dataset colour assignments."""
    js = _busco_chart_js(full_report)
    assert js.count("backgroundColor") == 4


def test_busco_chart_js_sparse_report_no_crash(sparse_report: GenomeReport) -> None:
    """_busco_chart_js does not crash when BUSCO data is None."""
    js = _busco_chart_js(sparse_report)
    assert "new Chart" in js
    assert "0.00" in js


# ---------------------------------------------------------------------------
# render_html — file creation
# ---------------------------------------------------------------------------


def test_render_html_creates_file(tmp_path: Path, full_report: GenomeReport) -> None:
    """render_html writes an HTML file to the output directory."""
    path = render_html(full_report, tmp_path)
    assert path.exists()


def test_render_html_filename_format(tmp_path: Path, full_report: GenomeReport) -> None:
    """The output filename replaces dots with underscores and ends in _report.html."""
    path = render_html(full_report, tmp_path)
    assert path.name == "GCA_000001515_5_report.html"


def test_render_html_returns_path_object(
    tmp_path: Path, full_report: GenomeReport
) -> None:
    """render_html returns a Path instance pointing to the written file."""
    path = render_html(full_report, tmp_path)
    assert isinstance(path, Path)


# ---------------------------------------------------------------------------
# render_html — HTML structure
# ---------------------------------------------------------------------------


def test_render_html_contains_species_name(
    tmp_path: Path, full_report: GenomeReport
) -> None:
    """The HTML report contains the scientific name in the header."""
    render_html(full_report, tmp_path)
    content = (tmp_path / "GCA_000001515_5_report.html").read_text(encoding="utf-8")
    assert "Pan troglodytes" in content


def test_render_html_contains_gca(tmp_path: Path, full_report: GenomeReport) -> None:
    """The HTML report contains the GCA accession string."""
    render_html(full_report, tmp_path)
    content = (tmp_path / "GCA_000001515_5_report.html").read_text(encoding="utf-8")
    assert "GCA_000001515.5" in content


def test_render_html_contains_quality_badge(
    tmp_path: Path, full_report: GenomeReport
) -> None:
    """The HTML report contains the quality badge text."""
    render_html(full_report, tmp_path)
    content = (tmp_path / "GCA_000001515_5_report.html").read_text(encoding="utf-8")
    assert "Excellent Quality" in content


def test_render_html_contains_busco_pct(
    tmp_path: Path, full_report: GenomeReport
) -> None:
    """The HTML report displays the BUSCO percentage in the metric card."""
    render_html(full_report, tmp_path)
    content = (tmp_path / "GCA_000001515_5_report.html").read_text(encoding="utf-8")
    assert "97.7%" in content


def test_render_html_contains_coding_genes(
    tmp_path: Path, full_report: GenomeReport
) -> None:
    """The HTML report displays the coding gene count."""
    render_html(full_report, tmp_path)
    content = (tmp_path / "GCA_000001515_5_report.html").read_text(encoding="utf-8")
    assert "23534" in content


def test_render_html_contains_ftp_link(
    tmp_path: Path, full_report: GenomeReport
) -> None:
    """The HTML report renders the FTP URL as a clickable anchor tag."""
    render_html(full_report, tmp_path)
    content = (tmp_path / "GCA_000001515_5_report.html").read_text(encoding="utf-8")
    assert "ftp.ebi.ac.uk" in content
    assert "<a href=" in content


def test_render_html_contains_chart_js_cdn(
    tmp_path: Path, full_report: GenomeReport
) -> None:
    """The HTML report loads Chart.js from the CDN."""
    render_html(full_report, tmp_path)
    content = (tmp_path / "GCA_000001515_5_report.html").read_text(encoding="utf-8")
    assert "cdn.jsdelivr.net" in content
    assert "chart.js" in content


def test_render_html_contains_footer(tmp_path: Path, full_report: GenomeReport) -> None:
    """The HTML report contains the Ensembl Genebuild Metadata footer."""
    render_html(full_report, tmp_path)
    content = (tmp_path / "GCA_000001515_5_report.html").read_text(encoding="utf-8")
    assert "Ensembl Genebuild Metadata" in content


def test_render_html_contains_metrics_table_headers(
    tmp_path: Path, full_report: GenomeReport
) -> None:
    """The metrics table contains expected row label text."""
    render_html(full_report, tmp_path)
    content = (tmp_path / "GCA_000001515_5_report.html").read_text(encoding="utf-8")
    assert "GCA Accession" in content
    assert "Coding Genes" in content
    assert "Assembly BUSCO" in content


def test_render_html_contains_busco_canvas(
    tmp_path: Path, full_report: GenomeReport
) -> None:
    """The HTML report contains a canvas element and Chart.js instantiation."""
    render_html(full_report, tmp_path)
    content = (tmp_path / "GCA_000001515_5_report.html").read_text(encoding="utf-8")
    assert 'id="buscoChart"' in content
    assert "new Chart" in content


# ---------------------------------------------------------------------------
# render_html — sparse / None field handling
# ---------------------------------------------------------------------------


def test_render_html_sparse_no_crash(
    tmp_path: Path, sparse_report: GenomeReport
) -> None:
    """render_html completes without error for a report with all-None optional fields."""
    path = render_html(sparse_report, tmp_path)
    assert path.exists()


def test_render_html_sparse_na_fallback(
    tmp_path: Path, sparse_report: GenomeReport
) -> None:
    """None field values are rendered as N/A in the HTML output."""
    render_html(sparse_report, tmp_path)
    content = (tmp_path / "GCA_000000001_1_report.html").read_text(encoding="utf-8")
    assert "N/A" in content


def test_render_html_sparse_no_ftp_anchor(
    tmp_path: Path, sparse_report: GenomeReport
) -> None:
    """When FTP is None, no anchor tag is rendered in the FTP row."""
    render_html(sparse_report, tmp_path)
    content = (tmp_path / "GCA_000000001_1_report.html").read_text(encoding="utf-8")
    assert "ftp.ebi.ac.uk" not in content


def test_render_html_busco_pct_zero_when_none(
    tmp_path: Path, sparse_report: GenomeReport
) -> None:
    """When protein_busco_complete is None, the displayed percentage is 0.0%."""
    render_html(sparse_report, tmp_path)
    content = (tmp_path / "GCA_000000001_1_report.html").read_text(encoding="utf-8")
    assert "0.0%" in content


def test_render_html_is_valid_html_skeleton(
    tmp_path: Path, full_report: GenomeReport
) -> None:
    """The generated file is a complete HTML document with all required tags."""
    render_html(full_report, tmp_path)
    content = (tmp_path / "GCA_000001515_5_report.html").read_text(encoding="utf-8")
    assert content.startswith("<!DOCTYPE html>")
    assert "</html>" in content
    assert "<body>" in content
    assert "</body>" in content
