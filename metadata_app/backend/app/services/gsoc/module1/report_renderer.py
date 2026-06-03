"""
report_renderer.py

Module 1: Per-genome report output generator.
Produces CSV, TXT, and PNG plot outputs from a GenomeReport.
"""

import csv
import logging
from dataclasses import asdict
from pathlib import Path
from typing import Optional

import matplotlib  # pylint: disable=wrong-import-order

matplotlib.use("Agg")  # Must be set before importing pyplot
import matplotlib.patches as mpatches  # noqa: E402  # pylint: disable=wrong-import-position
import matplotlib.pyplot as plt  # noqa: E402  # pylint: disable=wrong-import-position

from metadata_app.backend.app.services.gsoc.module1.busco_utils import (  # pylint: disable=wrong-import-position,import-error
    parse_busco_string,
)
from metadata_app.backend.app.services.gsoc.module1.genome_report import (  # pylint: disable=wrong-import-position,import-error
    GenomeReport,
)


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


def render_txt(report: GenomeReport, output_dir: Path) -> Path:
    """Write a human-readable text summary."""
    txt_path = output_dir / f"{report.gca.replace('.', '_')}_summary.txt"
    lines = [
        "=" * 60,
        "GENOME ANNOTATION QUALITY REPORT",
        "=" * 60,
        f"GCA Accession   : {report.gca}",
        f"Species         : {report.scientific_name}",
        f"Common Name     : {report.common_name or 'N/A'}",
        f"Taxon ID        : {report.lowest_taxon_id or 'N/A'}",
        f"Clade           : {report.internal_clade or 'N/A'}",
        "",
        "--- Annotation Info ---",
        f"Status          : {report.gb_status or 'N/A'}",
        f"Method          : {report.annotation_method or 'N/A'}",
        f"Genebuilder     : {report.genebuilder or 'N/A'}",
        f"Release Date    : {report.release_date or 'N/A'}",
        f"Latest Annotated: {report.latest_annotated or 'N/A'}",
        "",
        "--- Protein BUSCO ---",
        f"Score           : {report.protein_busco_raw or 'N/A'}",
        f"Complete        : {report.protein_busco_complete or 'N/A'}%",
        f"Quality         : {report.protein_busco_quality}",
        f"Lineage         : {report.protein_busco_lineage or 'N/A'}",
        "",
        "--- Assembly BUSCO ---",
        f"Score           : {report.assembly_busco_raw or 'N/A'}",
        f"Complete        : {report.assembly_busco_complete or 'N/A'}%",
        f"Lineage         : {report.assembly_busco_lineage or 'N/A'}",
        "",
        "--- Gene Statistics ---",
        f"Coding Genes    : {report.coding_genes or 'N/A'}",
        "",
        "--- FTP ---",
        f"Link            : {report.ftp or 'N/A'}",
        "=" * 60,
    ]
    with open(txt_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    logging.info("Written TXT to %s", txt_path)
    return txt_path


def _busco_color(complete: Optional[float]) -> str:
    """Return a hex color based on BUSCO completeness percentage."""
    if complete is None:
        return "#cccccc"
    if complete >= 95:
        return "#2ecc71"
    if complete >= 85:
        return "#f39c12"
    if complete >= 70:
        return "#e67e22"
    return "#e74c3c"


def _draw_busco_bars(ax: plt.Axes, raw: Optional[str], label: str) -> None:
    """Draw a stacked bar chart of BUSCO components onto the given axes."""
    parsed = parse_busco_string(raw or "")
    single = float(parsed.get("single_copy") or 0)
    duplicated = float(parsed.get("duplicated") or 0)
    fragmented = float(parsed.get("fragmented") or 0)
    missing = float(parsed.get("missing") or 0)
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
            (
                f"{report.protein_busco_complete}%"
                if report.protein_busco_complete
                else "N/A"
            ),
            _busco_color(report.protein_busco_complete),
        ),
        (
            "Assembly BUSCO",
            (
                f"{report.assembly_busco_complete}%"
                if report.assembly_busco_complete
                else "N/A"
            ),
            _busco_color(report.assembly_busco_complete),
        ),
        (
            "Quality Label",
            report.protein_busco_quality,
            _busco_color(report.protein_busco_complete),
        ),
        (
            "Coding Genes",
            str(report.coding_genes) if report.coding_genes else "N/A",
            "#3498db",
        ),
        ("Annotation Method", report.annotation_method or "N/A", "#9b59b6"),
        ("Status", report.gb_status or "N/A", "#1abc9c"),
        ("Clade", report.internal_clade or "N/A", "#e67e22"),
        ("Latest Annotated", report.latest_annotated or "N/A", "#2c3e50"),
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
