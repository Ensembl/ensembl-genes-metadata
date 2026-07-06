"""
comparison_renderer.py

Module 1: Multi-annotation comparison view.

For GCAs that have more than one genebuild_status row (e.g. a live
full_genebuild alongside an abandoned helixer attempt), this module
renders a side-by-side HTML table so reviewers can compare all
annotations at a glance.

Archive rows are excluded entirely per metrics_flagging.is_excluded_status().
All other rows (live, abandoned, completed, check_busco, etc.) are shown.

This module does NOT query the database directly - it receives the
anno_wide DataFrame produced by db_loader.load_anno_wide().
"""

import html
import logging
import math
import re
from datetime import date
from pathlib import Path
from typing import List

import pandas as pd

from metadata_app.backend.app.services.gsoc.module1.busco_utils import (  # pylint: disable=import-error
    busco_quality_label,
)
from metadata_app.backend.app.services.gsoc.module1.metrics_flagging import (  # pylint: disable=import-error
    is_excluded_status,
)

logger = logging.getLogger(__name__)

# Columns shown in the comparison table, in display order.
# Each entry is (column_name_in_dataframe, display_label).
_COMPARISON_COLUMNS = [
    ("gb_status", "Status"),
    ("annotation_method", "Method"),
    ("genebuilder", "Genebuilder"),
    ("annotated_version", "Version"),
    ("release_date", "Release Date"),
    ("date_status_update", "Last Status Update"),
    ("protein_busco", "Protein BUSCO"),
    ("protein_busco_lineage", "Protein Lineage"),
    ("assembly_busco", "Assembly BUSCO"),
    ("assembly_busco_lineage", "Assembly Lineage"),
    ("coding_genes", "Coding Genes"),
    ("annotation_source", "Source"),
    ("bioproject_id", "Bioproject"),
]

# Badge colours per gb_status for the column headers.
_STATUS_BADGE = {
    "live": ("color:#166534;background:#dcfce7", "Live"),
    "pre_released": ("color:#1e40af;background:#dbeafe", "Pre-released"),
    "handed_over": ("color:#5b21b6;background:#ede9fe", "Handed Over"),
    "coming_soon": ("color:#0f766e;background:#ccfbf1", "Coming Soon"),
    "completed": ("color:#92400e;background:#fef3c7", "Completed"),
    "check_busco": ("color:#9a3412;background:#ffedd5", "Check BUSCO"),
    "abandoned": ("color:#374151;background:#f3f4f6", "Abandoned"),
    "in_progress": ("color:#1d4ed8;background:#eff6ff", "In Progress"),
}
_DEFAULT_BADGE = ("color:#374151;background:#f3f4f6", "Unknown")

_CSS = """
    *, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
    body {
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
      background: #f8fafc;
      color: #1e293b;
      line-height: 1.6;
    }
    .header {
      background: #1e293b;
      color: #f8fafc;
      padding: 2rem 2.5rem;
    }
    .header h1 { font-size: 1.4rem; font-weight: 700; }
    .header .subtitle {
      font-size: 0.95rem; color: #94a3b8;
      margin-top: 0.2rem; font-style: italic;
    }
    .header .gca {
      font-size: 0.88rem; color: #64748b;
      margin-top: 0.15rem; font-family: monospace;
    }
    .container { max-width: 1100px; margin: 2rem auto; padding: 0 1.5rem; }
    .section {
      background: #fff; border: 1px solid #e2e8f0;
      border-radius: 0.75rem; padding: 1.5rem;
      margin-bottom: 1.5rem;
      box-shadow: 0 1px 3px rgba(0,0,0,0.06); overflow-x: auto;
    }
    .section h2 {
      font-size: 1rem; font-weight: 700; color: #1e293b;
      margin-bottom: 1rem; padding-bottom: 0.5rem;
      border-bottom: 1px solid #e2e8f0;
    }
    .footer {
      text-align: center; font-size: 0.8rem;
      color: #94a3b8; padding: 2rem 0;
    }
"""

_TABLE_STYLE = "width:100%;border-collapse:collapse;font-size:0.85rem;"


def _status_badge_style(gb_status: str) -> str:
    """Return inline CSS style string for a gb_status badge."""
    style, _ = _STATUS_BADGE.get(gb_status or "", _DEFAULT_BADGE)
    return style


def _status_badge_label(gb_status: str) -> str:
    """Return display label for a gb_status badge."""
    _, label = _STATUS_BADGE.get(gb_status or "", _DEFAULT_BADGE)
    return label


def _fmt(val: object) -> str:
    """Format a value for HTML display, replacing None/NaN/empty with N/A."""
    if val is None:
        return "N/A"
    if isinstance(val, float) and math.isnan(val):
        return "N/A"
    s = str(val).strip()
    return s if s and s.lower() != "none" else "N/A"


def _busco_cell_style(busco_raw: object, col: str) -> str:
    """Return a background colour hint for a BUSCO cell based on quality."""
    if "lineage" in col or "version" in col or "busco" not in col:
        return ""
    match = re.search(
        r"C:(\d+\.?\d*)%", str(busco_raw) if busco_raw is not None else ""
    )
    if not match:
        return ""
    label = busco_quality_label(float(match.group(1)))
    colors = {
        "Excellent": "background:#dcfce7",
        "Good": "background:#bbf7d0",
        "Moderate": "background:#fef3c7",
        "Poor": "background:#ffedd5",
    }
    return colors.get(label, "")


def _filter_rows(anno_wide: pd.DataFrame, gca: str) -> pd.DataFrame:
    """
    Return all non-archived annotation rows for a given GCA.

    Rows with gb_status == "archive" are excluded entirely per
    is_excluded_status(). All other statuses are included so the
    comparison view shows the full annotation history.
    """
    rows = anno_wide[anno_wide["gca"] == gca].copy()
    if rows.empty:
        return rows
    mask = rows["gb_status"].apply(
        lambda s: not is_excluded_status(str(s) if s is not None else "")
    )
    return rows[mask].reset_index(drop=True)


def _build_comparison_rows(rows: pd.DataFrame) -> List[tuple]:
    """
    Build a list of (display_label, col_name, [cell_value, ...]) tuples.

    Each tuple represents one metric row; the list of cell values has one
    entry per annotation column (i.e. one per non-archived genebuild_status
    row for this GCA).
    """
    result = []
    for col, label in _COMPARISON_COLUMNS:
        if col in rows.columns:
            values = [rows.iloc[i].get(col) for i in range(len(rows))]
        else:
            values = [None] * len(rows)
        result.append((label, col, values))
    return result


def _render_header_row(rows: pd.DataFrame) -> str:
    """Build the HTML table header row with one column per annotation."""
    th_style = (
        "padding:0.75rem;background:#f8fafc;text-align:center;"
        "border-bottom:2px solid #e2e8f0;min-width:160px;"
    )
    badge_wrap = (
        "display:inline-block;padding:0.2rem 0.7rem;"
        "border-radius:999px;font-size:0.78rem;font-weight:600;"
    )
    label_th = (
        "<th style='padding:0.75rem;background:#f8fafc;"
        "border-bottom:2px solid #e2e8f0;color:#475569;"
        "font-size:0.8rem;text-transform:uppercase;"
        "letter-spacing:0.05em'>Metric</th>"
    )
    cells = [label_th]
    for i in range(len(rows)):
        row = rows.iloc[i]
        status = str(row.get("gb_status") or "")
        method = _fmt(row.get("annotation_method"))
        version = _fmt(row.get("annotated_version"))
        badge_style = _status_badge_style(status)
        badge_label = _status_badge_label(status)
        cells.append(
            f"<th style='{th_style}'>"
            f"<span style='{badge_wrap}{html.escape(badge_style)}'>"
            f"{html.escape(badge_label)}</span>"
            f"<div style='font-size:0.82rem;color:#475569;margin-top:0.35rem'>"
            f"{html.escape(method)}</div>"
            f"<div style='font-size:0.78rem;color:#94a3b8'>"
            f"v{html.escape(version)}</div>"
            f"</th>"
        )
    return "<tr>" + "".join(cells) + "</tr>"


def _render_data_rows(comparison_rows: List[tuple]) -> str:
    """Build the HTML table data rows from pre-computed comparison tuples."""
    label_td_style = (
        "padding:0.55rem 0.75rem;font-weight:600;"
        "color:#475569;font-size:0.88rem;white-space:nowrap"
    )
    value_td_style = (
        "padding:0.55rem 0.75rem;font-size:0.88rem;font-family:monospace;color:#1e293b"
    )
    html_rows = []
    for idx, (label, col, values) in enumerate(comparison_rows):
        bg = "#ffffff" if idx % 2 == 0 else "#f8fafc"
        label_cell = f"<td style='{label_td_style}'>{html.escape(label)}</td>"
        value_cells = "".join(
            f"<td style='{value_td_style};{_busco_cell_style(v, col)}'>"
            f"{html.escape(_fmt(v))}</td>"
            for v in values
        )
        html_rows.append(f"<tr style='background:{bg}'>{label_cell}{value_cells}</tr>")
    return "".join(html_rows)


def render_comparison_html(
    gca: str,
    anno_wide: pd.DataFrame,
    output_dir: Path,
) -> Path:
    """
    Render a side-by-side HTML comparison table for all annotations of a GCA.

    Args:
        gca: GCA accession string e.g. "GCA_963455335.1"
        anno_wide: Full DataFrame returned by db_loader.load_anno_wide().
        output_dir: Directory to write the HTML file into.

    Returns:
        Path to the generated HTML file.

    Raises:
        ValueError: If the GCA is not found in anno_wide, or if all rows
            are excluded (archive only).
    """
    rows = _filter_rows(anno_wide, gca)
    if rows.empty:
        raise ValueError(
            f"GCA '{gca}' not found in anno_wide or all rows are excluded."
        )

    logger.info(
        "Rendering comparison view for GCA %s (%d annotation rows)", gca, len(rows)
    )

    scientific_name = _fmt(rows.iloc[0].get("scientific_name"))
    common_name = _fmt(rows.iloc[0].get("common_name"))
    generated = date.today().isoformat()
    subtitle = scientific_name
    if common_name != "N/A":
        subtitle += f" ({common_name})"

    header_row = _render_header_row(rows)
    comparison_rows = _build_comparison_rows(rows)
    data_rows = _render_data_rows(comparison_rows)
    table_html = (
        f"<table style='{_TABLE_STYLE}'>"
        f"<thead>{header_row}</thead>"
        f"<tbody>{data_rows}</tbody>"
        f"</table>"
    )

    html_doc = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8"/>
  <meta name="viewport" content="width=device-width, initial-scale=1.0"/>
  <title>Annotation Comparison — {html.escape(gca)}</title>
  <style>{_CSS}</style>
</head>
<body>
<div class="header">
  <h1>Annotation Comparison</h1>
  <div class="subtitle">{html.escape(subtitle)}</div>
  <div class="gca">{html.escape(gca)} &mdash; {len(rows)} annotation(s)</div>
</div>
<div class="container">
  <div class="section">
    <h2>Side-by-Side Comparison</h2>
    {table_html}
  </div>
  <div class="footer">
    Generated by Ensembl Genebuild Metadata &mdash; {generated}
  </div>
</div>
</body>
</html>
"""

    html_path = output_dir / f"{gca.replace('.', '_')}_comparison.html"
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(html_doc)

    logger.info("Written comparison HTML to %s", html_path)
    return html_path
