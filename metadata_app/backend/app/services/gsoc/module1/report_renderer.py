"""
report_renderer.py

Module 1: Per-genome report output generator.
Produces CSV, TXT, and PNG plot outputs from a GenomeReport.
"""

import csv
import logging
from dataclasses import asdict
from pathlib import Path
from typing import Dict, Optional, Union

import matplotlib  # pylint: disable=wrong-import-order

matplotlib.use("Agg")  # Must be set before importing pyplot
import matplotlib.patches as mpatches  # noqa: E402  # pylint: disable=wrong-import-position
import matplotlib.pyplot as plt  # noqa: E402  # pylint: disable=wrong-import-position

from metadata_app.backend.app.services.gsoc.module1.busco_utils import (  # pylint: disable=wrong-import-position,import-error
    QUALITY_THRESHOLDS,
    parse_busco_string,
)
from metadata_app.backend.app.services.gsoc.module1.genome_report import (  # pylint: disable=wrong-import-position,import-error
    GenomeReport,
)

# Grey used for "no data" / unknown across all BUSCO-related visuals.
_NO_DATA_COLOR = "#cccccc"

# Colour for the worst quality band, below the lowest QUALITY_THRESHOLDS entry.
_POOR_COLOR = "#e74c3c"

# Explicit colour-per-band mapping, in the same order as QUALITY_THRESHOLDS
# (highest threshold first). Keeping this as its own small table (rather
# than re-deriving colours mathematically) makes the colour choices easy
# to tweak without touching the threshold logic itself.
_BAND_COLORS = ("#2ecc71", "#f39c12", "#e67e22")


def create_output_directory(gca: str, base_output_dir: str = "outputs") -> Path:
    """Create output directory for this genome report."""
    safe_gca = gca.replace(".", "_")
    output_path = Path(base_output_dir) / safe_gca
    output_path.mkdir(parents=True, exist_ok=True)
    logging.info("Created output directory: %s", output_path)
    return output_path


def render_csv(report: GenomeReport, output_dir: Path) -> Path:
    """Write all GenomeReport metrics to a CSV file."""
    csv_path = output_dir / f"{report.gca.replace('.', '_')}_metrics.csv"
    report_dict = asdict(report)
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "value"])
        for key, value in report_dict.items():
            writer.writerow([key, value if value is not None else ""])
    logging.info("Written CSV to %s", csv_path)
    return csv_path


def _fmt_value(value: Optional[object]) -> str:
    """
    Format an optional field for display, using an explicit None check.

    Using `value is not None` (rather than truthiness) matters here because
    a legitimate value of 0, 0.0, or "" must still be displayed as-is, not
    silently replaced with "N/A". Falling back to truthiness would, for
    example, show "N/A" for a genuinely 0% BUSCO score.
    """
    return str(value) if value is not None else "N/A"


def _fmt_busco_pct(value: Optional[float]) -> str:
    """
    Format a BUSCO completeness percentage consistently to 1 decimal place.

    Uses an explicit None check (not truthiness) so a real 0.0% value is
    still displayed as "0.0%" rather than falling back to "N/A".
    """
    if value is None:
        return "N/A"
    return f"{value:.1f}%"


def render_txt(report: GenomeReport, output_dir: Path) -> Path:
    """Write a human-readable text summary."""
    txt_path = output_dir / f"{report.gca.replace('.', '_')}_summary.txt"
    lines = [
        "=" * 60,
        "GENOME ANNOTATION QUALITY REPORT",
        "=" * 60,
        f"GCA Accession   : {report.gca}",
        f"Species         : {report.scientific_name}",
        f"Common Name     : {_fmt_value(report.common_name)}",
        f"Taxon ID        : {_fmt_value(report.lowest_taxon_id)}",
        f"Clade           : {_fmt_value(report.internal_clade)}",
        "",
        "--- Annotation Info ---",
        f"Status          : {_fmt_value(report.gb_status)}",
        f"Method          : {_fmt_value(report.annotation_method)}",
        f"Genebuilder     : {_fmt_value(report.genebuilder)}",
        f"Release Date    : {_fmt_value(report.release_date)}",
        f"Latest Annotated: {_fmt_value(report.latest_annotated)}",
        "",
        "--- Protein BUSCO ---",
        f"Score           : {_fmt_value(report.protein_busco_raw)}",
        f"Complete        : {_fmt_busco_pct(report.protein_busco_complete)}",
        f"Quality         : {report.protein_busco_quality}",
        f"Lineage         : {_fmt_value(report.protein_busco_lineage)}",
        "",
        "--- Assembly BUSCO ---",
        f"Score           : {_fmt_value(report.assembly_busco_raw)}",
        f"Complete        : {_fmt_busco_pct(report.assembly_busco_complete)}",
        f"Lineage         : {_fmt_value(report.assembly_busco_lineage)}",
        "",
        "--- Gene Statistics ---",
        f"Coding Genes    : {_fmt_value(report.coding_genes)}",
        "",
        "--- FTP ---",
        f"Link            : {_fmt_value(report.ftp)}",
        "=" * 60,
    ]
    with open(txt_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    logging.info("Written TXT to %s", txt_path)
    return txt_path


def _busco_color(complete: Optional[float]) -> str:
    """
    Return a hex color based on BUSCO completeness percentage.

    Sources thresholds from busco_utils.QUALITY_THRESHOLDS (the single
    canonical source of truth for BUSCO quality bands) rather than
    redefining the 95/85/70 cutoffs locally, so this stays in sync with
    busco_quality_label() and cannot silently drift out of step with it.
    """
    if complete is None:
        return _NO_DATA_COLOR
    for (_, threshold), color in zip(QUALITY_THRESHOLDS, _BAND_COLORS):
        if complete >= threshold:
            return color
    return _POOR_COLOR


def _numeric_field(
    parsed: Dict[str, Optional[Union[float, int, Dict[str, float]]]], key: str
) -> float:
    """
    Safely extract a numeric (float or int) value from a parsed BUSCO dict.

    parse_busco_string() returns Optional[Union[float, int, Dict[str, float]]]
    per field because the same return type is shared with the "extra" key,
    which holds a dict. Standard fields (single_copy, duplicated, fragmented,
    missing, etc.) never actually hold a dict in practice, but the type
    system can't express that distinction, so we narrow explicitly here
    with isinstance() rather than silencing mypy with type:ignore. Returns
    0.0 for None, missing keys, or any unexpected non-numeric value.
    """
    value = parsed.get(key)
    if isinstance(value, (float, int)):
        return float(value)
    return 0.0


def _draw_busco_bars(ax: plt.Axes, raw: Optional[str], label: str) -> None:
    """Draw a stacked bar chart of BUSCO components onto the given axes."""
    parsed = parse_busco_string(raw or "")
    single = _numeric_field(parsed, "single_copy")
    duplicated = _numeric_field(parsed, "duplicated")
    fragmented = _numeric_field(parsed, "fragmented")
    missing = _numeric_field(parsed, "missing")
    total = single + duplicated + fragmented + missing
    if total == 0:
        ax.text(0.5, 0.5, "No BUSCO data", ha="center", va="center")
        ax.set_title(label)
        return
    components = [
        (single, "#2ecc71", "Single"),
        (duplicated, "#27ae60", "Duplicated"),
        (fragmented, "#f39c12", "Fragmented"),
        (missing, "#e74c3c", "Missing"),
    ]
    bottom: float = 0.0
    for val, color, bar_label in components:
        ax.bar(0, val, bottom=bottom, color=color, width=0.5, label=bar_label)
        if val > 1:
            ax.text(
                0,
                bottom + val / 2,
                f"{val:.1f}%",
                ha="center",
                va="center",
                fontsize=9,
                color="white",
                fontweight="bold",
            )
        bottom += val
    ax.set_title(label, fontsize=11)
    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(0, 105)
    ax.set_ylabel("Percentage (%)")
    ax.set_xticks([])
    ax.legend(loc="upper right", fontsize=8)


def plot_busco_bar(report: GenomeReport, output_dir: Path) -> Path:
    """Plot stacked bar chart of BUSCO composition for protein and assembly."""
    plot_path = output_dir / f"{report.gca.replace('.', '_')}_busco.png"
    _, axes = plt.subplots(1, 2, figsize=(12, 5))
    plt.suptitle(
        f"BUSCO Scores — {report.scientific_name} ({report.gca})",
        fontsize=13,
        fontweight="bold",
    )
    _draw_busco_bars(axes[0], report.protein_busco_raw, "Protein BUSCO")
    _draw_busco_bars(axes[1], report.assembly_busco_raw, "Assembly BUSCO")
    plt.tight_layout()
    plt.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close()
    logging.info("Saved BUSCO plot to %s", plot_path)
    return plot_path


def plot_quality_summary(report: GenomeReport, output_dir: Path) -> Path:
    """Plot a quality summary card showing key metrics at a glance."""
    plot_path = output_dir / f"{report.gca.replace('.', '_')}_quality_summary.png"
    _, ax = plt.subplots(figsize=(8, 5))
    ax.axis("off")
    ax.text(
        0.5,
        0.95,
        f"{report.scientific_name}\n{report.gca}",
        transform=ax.transAxes,
        fontsize=14,
        fontweight="bold",
        ha="center",
        va="top",
    )
    metrics = [
        (
            "Protein BUSCO",
            _fmt_busco_pct(report.protein_busco_complete),
            _busco_color(report.protein_busco_complete),
        ),
        (
            "Assembly BUSCO",
            _fmt_busco_pct(report.assembly_busco_complete),
            _busco_color(report.assembly_busco_complete),
        ),
        (
            "Quality Label",
            report.protein_busco_quality,
            _busco_color(report.protein_busco_complete),
        ),
        (
            "Coding Genes",
            _fmt_value(report.coding_genes),
            "#3498db",
        ),
        ("Annotation Method", _fmt_value(report.annotation_method), "#9b59b6"),
        ("Status", _fmt_value(report.gb_status), "#1abc9c"),
        ("Clade", _fmt_value(report.internal_clade), "#e67e22"),
        ("Latest Annotated", _fmt_value(report.latest_annotated), "#2c3e50"),
    ]
    y = 0.82
    for metric_name, value, color in metrics:
        ax.add_patch(
            mpatches.FancyBboxPatch(
                (0.02, y - 0.04),
                0.96,
                0.08,
                boxstyle="round,pad=0.01",
                facecolor=color,
                alpha=0.15,
                transform=ax.transAxes,
                clip_on=False,
            )
        )
        ax.text(
            0.05,
            y,
            metric_name,
            transform=ax.transAxes,
            fontsize=10,
            va="center",
            color="#2c3e50",
        )
        ax.text(
            0.95,
            y,
            value,
            transform=ax.transAxes,
            fontsize=10,
            va="center",
            ha="right",
            fontweight="bold",
            color=color,
        )
        y -= 0.10
    plt.tight_layout()
    plt.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close()
    logging.info("Saved quality summary plot to %s", plot_path)
    return plot_path


def render_report(report: GenomeReport, base_output_dir: str = "outputs") -> dict:
    """Main entry point. Generate all outputs for a GenomeReport."""
    logging.info("Rendering report for %s", report.gca)
    output_dir = create_output_directory(report.gca, base_output_dir)
    csv_path = render_csv(report, output_dir)
    txt_path = render_txt(report, output_dir)
    busco_plot = plot_busco_bar(report, output_dir)
    summary_plot = plot_quality_summary(report, output_dir)
    outputs = {
        "csv": str(csv_path),
        "txt": str(txt_path),
        "busco_plot": str(busco_plot),
        "summary_plot": str(summary_plot),
        "output_dir": str(output_dir),
    }
    logging.info("Report complete for %s", report.gca)
    return outputs
