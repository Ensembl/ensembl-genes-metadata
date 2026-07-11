"""
test_report_renderer.py

Unit tests for report_renderer.py.
All tests use synthetic fixture data — no database connection required.
"""

# pylint: disable=redefined-outer-name,import-outside-toplevel

import csv
from pathlib import Path

import pytest

from metadata_app.backend.app.services.gsoc.module1.genome_report import GenomeReport
from metadata_app.backend.app.services.gsoc.module1.report_renderer import (
    _busco_color,
    _draw_busco_bars,
    create_output_directory,
    plot_busco_bar,
    plot_quality_summary,
    render_csv,
    render_report,
    render_txt,
)

# Fixtures


@pytest.fixture()
def full_report() -> GenomeReport:
    """A fully-populated GenomeReport for a mock human genome."""
    return GenomeReport(
        gca="GCA_000001405.29",
        scientific_name="Homo sapiens",
        common_name="human",
        lowest_taxon_id=9606,
        internal_clade="Primates",
        annotation_method="full_genebuild",
        genebuilder="genebuild_team",
        gb_status="live",
        annotation_source="ensembl",
        bioproject_id="PRJNA31257",
        associated_project="Human Genome Project",
        release_date="2023-01-15",
        date_status_update="2023-01-15",
        last_genebuild_update="2023-01-10",
        protein_busco_raw="C:94.3%[S:91.2%,D:3.1%],F:2.1%,M:3.6%,n:255",
        protein_busco_lineage="primates_odb10",
        protein_busco_version="5.4.3",
        protein_busco_complete=94.3,
        protein_busco_quality="High",
        protein_busco_extra={},
        assembly_busco_raw="C:98.1%[S:97.0%,D:1.1%],F:0.9%,M:1.0%,n:255",
        assembly_busco_lineage="primates_odb10",
        assembly_busco_version="5.4.3",
        assembly_busco_complete=98.1,
        assembly_busco_extra={},
        coding_genes=20442,
        total_transcripts=61225,
        transcripts_per_gene=1.83,
        average_cds_length=1603.0,
        average_coding_intron_length=5560.0,
        single_exon_coding_genes=2190,
        longest_coding_gene_length=2513825,
        average_coding_exon_length=157.0,
        nc_non_coding_genes=9710,
        latest_annotated="Yes",
        annotated_version=29.0,
        assembly_version=29.0,
        ftp="https://ftp.ebi.ac.uk/pub/ensemblorganisms/Homo_sapiens/GCA_000001405.29/",
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
        protein_busco_quality="Unknown",
        protein_busco_extra={},
        assembly_busco_raw=None,
        assembly_busco_lineage=None,
        assembly_busco_version=None,
        assembly_busco_complete=None,
        assembly_busco_extra={},
        coding_genes=None,
        total_transcripts=None,
        transcripts_per_gene=None,
        average_cds_length=None,
        average_coding_intron_length=None,
        single_exon_coding_genes=None,
        longest_coding_gene_length=None,
        average_coding_exon_length=None,
        nc_non_coding_genes=None,
        latest_annotated=None,
        annotated_version=None,
        assembly_version=None,
        ftp=None,
    )


# create_output_directory


def test_create_output_directory_creates_path(tmp_path: Path) -> None:
    """Directory is created and GCA dots are replaced with underscores."""
    result = create_output_directory("GCA_000001405.29", str(tmp_path))
    assert result.exists()
    assert result.is_dir()
    assert result.name == "GCA_000001405_29"


def test_create_output_directory_idempotent(tmp_path: Path) -> None:
    """Calling twice with the same GCA does not raise."""
    create_output_directory("GCA_000001405.29", str(tmp_path))
    create_output_directory("GCA_000001405.29", str(tmp_path))


# render_csv


def test_render_csv_creates_file(tmp_path: Path, full_report: GenomeReport) -> None:
    """CSV file is created in the expected location."""
    csv_path = render_csv(full_report, tmp_path)
    assert csv_path.exists()
    assert csv_path.suffix == ".csv"


def test_render_csv_has_header(tmp_path: Path, full_report: GenomeReport) -> None:
    """CSV first row is the metric/value header."""
    csv_path = render_csv(full_report, tmp_path)
    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
    assert header == ["metric", "value"]


def test_render_csv_contains_gca(tmp_path: Path, full_report: GenomeReport) -> None:
    """CSV contains a row with the GCA accession."""
    csv_path = render_csv(full_report, tmp_path)
    with open(csv_path, encoding="utf-8") as f:
        content = f.read()
    assert "GCA_000001405.29" in content


def test_render_csv_none_values_become_empty(
    tmp_path: Path, sparse_report: GenomeReport
) -> None:
    """None values are written as empty strings, not the word None."""
    csv_path = render_csv(sparse_report, tmp_path)
    with open(csv_path, encoding="utf-8") as f:
        content = f.read()
    assert "None" not in content


# render_txt


def test_render_txt_creates_file(tmp_path: Path, full_report: GenomeReport) -> None:
    """TXT file is created in the expected location."""
    txt_path = render_txt(full_report, tmp_path)
    assert txt_path.exists()
    assert txt_path.suffix == ".txt"


def test_render_txt_contains_species(tmp_path: Path, full_report: GenomeReport) -> None:
    """TXT summary contains the scientific name."""
    txt_path = render_txt(full_report, tmp_path)
    assert "Homo sapiens" in txt_path.read_text(encoding="utf-8")


def test_render_txt_contains_gca(tmp_path: Path, full_report: GenomeReport) -> None:
    """TXT summary contains the GCA accession."""
    txt_path = render_txt(full_report, tmp_path)
    assert "GCA_000001405.29" in txt_path.read_text(encoding="utf-8")


def test_render_txt_sparse_no_crash(
    tmp_path: Path, sparse_report: GenomeReport
) -> None:
    """render_txt does not crash on a report with all-None optional fields."""
    txt_path = render_txt(sparse_report, tmp_path)
    assert txt_path.exists()


# _busco_color


def test_busco_color_none_returns_grey() -> None:
    """None input returns the no-data grey colour."""
    assert _busco_color(None) == "#cccccc"


def test_busco_color_high() -> None:
    """A score of 95 returns the green colour."""
    assert _busco_color(95.0) == "#2ecc71"


def test_busco_color_medium_high() -> None:
    """A score of 90 returns the amber colour."""
    assert _busco_color(90.0) == "#f39c12"


def test_busco_color_medium_low() -> None:
    """A score of 75 returns the orange colour."""
    assert _busco_color(75.0) == "#e67e22"


def test_busco_color_low() -> None:
    """A score of 60 returns the red colour."""
    assert _busco_color(60.0) == "#e74c3c"


def test_busco_color_boundary_95() -> None:
    """Exactly 95 should return green."""
    assert _busco_color(95.0) == "#2ecc71"


def test_busco_color_boundary_85() -> None:
    """Exactly 85 should return amber."""
    assert _busco_color(85.0) == "#f39c12"


# plot_busco_bar


def test_plot_busco_bar_creates_file(tmp_path: Path, full_report: GenomeReport) -> None:
    """BUSCO bar plot PNG is created."""
    plot_path = plot_busco_bar(full_report, tmp_path)
    assert plot_path.exists()
    assert plot_path.suffix == ".png"


def test_plot_busco_bar_sparse_no_crash(
    tmp_path: Path, sparse_report: GenomeReport
) -> None:
    """BUSCO bar plot does not crash when BUSCO data is None."""
    plot_path = plot_busco_bar(sparse_report, tmp_path)
    assert plot_path.exists()


# plot_quality_summary


def test_plot_quality_summary_creates_file(
    tmp_path: Path, full_report: GenomeReport
) -> None:
    """Quality summary PNG is created."""
    plot_path = plot_quality_summary(full_report, tmp_path)
    assert plot_path.exists()
    assert plot_path.suffix == ".png"


def test_plot_quality_summary_sparse_no_crash(
    tmp_path: Path, sparse_report: GenomeReport
) -> None:
    """Quality summary does not crash when optional fields are None."""
    plot_path = plot_quality_summary(sparse_report, tmp_path)
    assert plot_path.exists()


# _draw_busco_bars (internal helper)


def test_draw_busco_bars_no_crash_empty_string() -> None:
    """_draw_busco_bars handles empty BUSCO string without crashing."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _, ax = plt.subplots()
    _draw_busco_bars(ax, "", "Test")
    plt.close()


def test_draw_busco_bars_no_crash_none() -> None:
    """_draw_busco_bars handles None BUSCO string without crashing."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _, ax = plt.subplots()
    _draw_busco_bars(ax, None, "Test")
    plt.close()


# render_report (integration)


def test_render_report_returns_all_keys(
    tmp_path: Path, full_report: GenomeReport
) -> None:
    """render_report returns a dict with all expected output keys."""
    outputs = render_report(full_report, base_output_dir=str(tmp_path))
    assert set(outputs.keys()) == {
        "csv",
        "txt",
        "busco_plot",
        "summary_plot",
        "output_dir",
    }


def test_render_report_all_files_exist(
    tmp_path: Path, full_report: GenomeReport
) -> None:
    """All files referenced in the render_report output actually exist."""
    outputs = render_report(full_report, base_output_dir=str(tmp_path))
    for key in ("csv", "txt", "busco_plot", "summary_plot"):
        assert Path(outputs[key]).exists(), f"Missing: {key}"


def test_render_report_sparse_no_crash(
    tmp_path: Path, sparse_report: GenomeReport
) -> None:
    """render_report completes without error for a sparse report."""
    outputs = render_report(sparse_report, base_output_dir=str(tmp_path))
    assert Path(outputs["csv"]).exists()
