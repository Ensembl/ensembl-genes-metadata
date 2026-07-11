"""
html_renderer.py

Module 1: HTML report generator.
Produces a self-contained HTML quality report from a GenomeReport.
No external dependencies at render time - Chart.js loaded from CDN.
"""

import html
import logging
from datetime import date
from pathlib import Path
from typing import Dict, Optional, Union

from metadata_app.backend.app.services.gsoc.module1.busco_utils import (  # pylint: disable=import-error
    QUALITY_POOR,
    QUALITY_THRESHOLDS,
    QUALITY_UNKNOWN,
    parse_busco_string,
)
from metadata_app.backend.app.services.gsoc.module1.genome_report import (  # pylint: disable=import-error
    GenomeReport,
)

logger = logging.getLogger(__name__)

# Colour pairs (text colour, background colour), one per QUALITY_THRESHOLDS
# band in order (highest threshold first), plus explicit entries for the
# Poor/Unknown bands which fall outside QUALITY_THRESHOLDS itself.
_BAND_BADGE_COLORS = (
    ("#166534", "#dcfce7"),  # Excellent
    ("#14532d", "#bbf7d0"),  # Good
    ("#92400e", "#fef3c7"),  # Moderate
)
_POOR_BADGE_COLOR = ("#9a3412", "#ffedd5")
_UNKNOWN_BADGE_COLOR = ("#374151", "#f3f4f6")

# Same band colours, but as flat hex strings for the progress-bar fill
# (no background pairing needed there).
_BAND_BAR_COLORS = ("#16a34a", "#d97706", "#ea580c")
_POOR_BAR_COLOR = "#dc2626"
_NO_DATA_BAR_COLOR = "#9ca3af"


def _quality_badge_style(quality: str) -> str:
    """
    Return CSS colour pair for a quality label.

    Sources its label set from busco_utils (QUALITY_THRESHOLDS band names,
    plus QUALITY_POOR/QUALITY_UNKNOWN) so this stays in sync with the
    canonical vocabulary rather than redefining its own label strings,
    which previously caused every label except "Excellent" to silently
    fall back to the default grey.
    """
    band_names = [name for name, _ in QUALITY_THRESHOLDS]
    mapping: Dict[str, tuple] = dict(zip(band_names, _BAND_BADGE_COLORS))
    mapping[QUALITY_POOR] = _POOR_BADGE_COLOR
    mapping[QUALITY_UNKNOWN] = _UNKNOWN_BADGE_COLOR
    color, bg = mapping.get(quality, _UNKNOWN_BADGE_COLOR)
    return f"color:{color};background:{bg}"


def _busco_bar_color(pct: Optional[float]) -> str:
    """
    Return a hex fill colour for the BUSCO progress bar.

    Sources thresholds from busco_utils.QUALITY_THRESHOLDS rather than
    redefining the 95/85/70 cutoffs locally, so this cannot silently
    drift out of sync with busco_quality_label() or report_renderer's
    _busco_color().
    """
    if pct is None:
        return _NO_DATA_BAR_COLOR
    for (_, threshold), color in zip(QUALITY_THRESHOLDS, _BAND_BAR_COLORS):
        if pct >= threshold:
            return color
    return _POOR_BAR_COLOR


def _fmt(val: object) -> str:
    """Format a value for display, replacing None/empty with N/A."""
    if val is None or val == "":
        return "N/A"
    return str(val)


def _numeric_field(
    parsed: Dict[str, Optional[Union[float, int, Dict[str, float]]]], key: str
) -> float:
    """
    Safely extract a numeric (float or int) value from a parsed BUSCO dict.

    parse_busco_string() returns Optional[Union[float, int, Dict[str, float]]]
    per field because the same return type is shared with the "extra" key,
    which holds a dict. Standard fields never actually hold a dict in
    practice, but the type system can't express that, so we narrow
    explicitly with isinstance() rather than silencing mypy with
    type:ignore. Returns 0.0 for None, missing keys, or any unexpected
    non-numeric value.
    """
    value = parsed.get(key)
    if isinstance(value, (float, int)):
        return float(value)
    return 0.0


def _safe_float(v: object) -> float:
    """
    Safely cast a value to float, returning 0.0 on failure.

    Explicitly excludes None and dict (the latter being the type
    parse_busco_string()'s "extra" field can hold) before attempting the
    cast, then relies on float()'s own str/int/float handling for
    everything else. The remaining type:ignore is narrow and deliberate:
    it covers only "trust that a non-None, non-dict object is float()-able
    at runtime", not the dict-confusion bug this replaced.
    """
    if v is None or isinstance(v, dict):
        return 0.0
    try:
        return float(v)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return 0.0


def _busco_dataset_row(
    parsed: Dict[str, Optional[Union[float, int, Dict[str, float]]]],
) -> tuple:
    """Return (single, duplicated, fragmented, missing) as floats for one BUSCO row."""
    return (
        _numeric_field(parsed, "single_copy"),
        _numeric_field(parsed, "duplicated"),
        _numeric_field(parsed, "fragmented"),
        _numeric_field(parsed, "missing"),
    )


def _busco_chart_js(report: GenomeReport) -> str:
    """
    Return an inline Chart.js script block for the protein + assembly
    BUSCO stacked bars.

    If neither protein nor assembly BUSCO data is available, the chart
    is still rendered but with an explicit "No BUSCO data available"
    message overlay instead of silently showing an empty/invisible chart.
    """
    protein_row = _busco_dataset_row(parse_busco_string(report.protein_busco_raw or ""))
    assembly_row = _busco_dataset_row(
        parse_busco_string(report.assembly_busco_raw or "")
    )

    if (sum(protein_row) + sum(assembly_row)) == 0:
        return """
    const ctx = document.getElementById('buscoChart').getContext('2d');
    ctx.canvas.parentElement.innerHTML =
        '<p style="text-align:center;color:#94a3b8;padding:2rem 0;">' +
        'No BUSCO data available</p>';
"""

    # (label, colour, index-into-protein_row/assembly_row) for each BUSCO
    # component, so the four datasets can be built in one loop instead of
    # four near-identical blocks of locals.
    component_specs = (
        ("Single", "#16a34a", 0),
        ("Duplicated", "#4ade80", 1),
        ("Fragmented", "#f59e0b", 2),
        ("Missing", "#ef4444", 3),
    )
    datasets = [
        (
            f"{{ label: '{label}', data: [{protein_row[idx]:.2f}, {assembly_row[idx]:.2f}],"
            f" backgroundColor: '{color}' }}"
        )
        for label, color, idx in component_specs
    ]
    ds_single, ds_dup, ds_frag, ds_missing = datasets

    return f"""
    const ctx = document.getElementById('buscoChart').getContext('2d');
    new Chart(ctx, {{
        type: 'bar',
        data: {{
            labels: ['Protein BUSCO', 'Assembly BUSCO'],
            datasets: [
                {ds_single},
                {ds_dup},
                {ds_frag},
                {ds_missing},
            ]
        }},
        options: {{
            indexAxis: 'y',
            responsive: true,
            plugins: {{
                legend: {{ position: 'bottom', labels: {{ font: {{ size: 12 }} }} }},
                title: {{ display: false }}
            }},
            scales: {{
                x: {{
                    stacked: true,
                    max: 100,
                    title: {{ display: true, text: 'Percentage (%)' }}
                }},
                y: {{ stacked: true }}
            }}
        }}
    }});
"""


def render_html(  # pylint: disable=too-many-locals
    report: GenomeReport, output_dir: Path
) -> Path:
    """
    Render a self-contained HTML quality report for a GenomeReport.

    Note: many local variables here are inherent to assembling one large
    HTML document from many independent report fields (badge style, two
    BUSCO percentages/colours, chart script, metrics rows, ftp html, four
    display strings, etc.) -- grouping them into a dict would reduce the
    local count but make the long template harder to read, since each
    name documents what it is at the point of use.

    Args:
        report: Populated GenomeReport dataclass.
        output_dir: Directory to write the HTML file into.

    Returns:
        Path to the generated HTML file.
    """
    html_path = output_dir / f"{report.gca.replace('.', '_')}_report.html"
    badge_style = _quality_badge_style(report.protein_busco_quality)
    protein_busco_pct = report.protein_busco_complete
    assembly_busco_pct = report.assembly_busco_complete
    protein_bar_color = _busco_bar_color(protein_busco_pct)
    assembly_bar_color = _busco_bar_color(assembly_busco_pct)
    chart_script = _busco_chart_js(report)
    generated = date.today().isoformat()

    metrics_rows = [
        ("GCA Accession", _fmt(report.gca)),
        ("Scientific Name", _fmt(report.scientific_name)),
        ("Common Name", _fmt(report.common_name)),
        ("Taxon ID", _fmt(report.lowest_taxon_id)),
        ("Clade", _fmt(report.internal_clade)),
        ("Annotation Status", _fmt(report.gb_status)),
        ("Annotation Method", _fmt(report.annotation_method)),
        ("Genebuilder", _fmt(report.genebuilder)),
        ("Annotation Source", _fmt(report.annotation_source)),
        ("Bioproject ID", _fmt(report.bioproject_id)),
        ("Release Date", _fmt(report.release_date)),
        ("Date Status Update", _fmt(report.date_status_update)),
        ("Last Genebuild Update", _fmt(report.last_genebuild_update)),
        ("Protein BUSCO", _fmt(report.protein_busco_raw)),
        ("Protein BUSCO Lineage", _fmt(report.protein_busco_lineage)),
        ("Protein BUSCO Version", _fmt(report.protein_busco_version)),
        ("Assembly BUSCO", _fmt(report.assembly_busco_raw)),
        ("Assembly BUSCO Lineage", _fmt(report.assembly_busco_lineage)),
        ("Assembly BUSCO Version", _fmt(report.assembly_busco_version)),
        ("Coding Genes", _fmt(report.coding_genes)),
        ("Total Transcripts", _fmt(report.total_transcripts)),
        ("Transcripts per Gene", _fmt(report.transcripts_per_gene)),
        ("Avg CDS Length (bp)", _fmt(report.average_cds_length)),
        ("Avg Coding Intron Length (bp)", _fmt(report.average_coding_intron_length)),
        ("Single Exon Coding Genes", _fmt(report.single_exon_coding_genes)),
        ("Longest Coding Gene (bp)", _fmt(report.longest_coding_gene_length)),
        ("Avg Coding Exon Length (bp)", _fmt(report.average_coding_exon_length)),
        ("Non-Coding Genes", _fmt(report.nc_non_coding_genes)),
        ("Latest Annotated", _fmt(report.latest_annotated)),
        ("Annotated Version", _fmt(report.annotated_version)),
        ("Assembly Version", _fmt(report.assembly_version)),
    ]

    table_rows_html = "\n".join(
        f"            <tr>\n"
        f'              <td class="metric-name">{html.escape(name)}</td>\n'
        f'              <td class="metric-value">{html.escape(value)}</td>\n'
        f"            </tr>"
        for name, value in metrics_rows
    )

    if report.ftp:
        safe_ftp = html.escape(report.ftp, quote=True)
        ftp_html = (
            f'<a href="{safe_ftp}" target="_blank" class="ftp-link">{safe_ftp}</a>'
        )
    else:
        ftp_html = "N/A"

    protein_bar_pct = (
        f"{min(protein_busco_pct, 100):.1f}%" if protein_busco_pct is not None else "0%"
    )
    assembly_bar_pct = (
        f"{min(assembly_busco_pct, 100):.1f}%"
        if assembly_busco_pct is not None
        else "0%"
    )
    protein_pct_display = (
        f"{protein_busco_pct:.1f}%" if protein_busco_pct is not None else "N/A"
    )
    assembly_pct_display = (
        f"{assembly_busco_pct:.1f}%" if assembly_busco_pct is not None else "N/A"
    )

    html_doc = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8"/>
  <meta name="viewport" content="width=device-width, initial-scale=1.0"/>
  <title>Genome Report — {html.escape(report.scientific_name)} ({html.escape(report.gca)})</title>
  <script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
  <style>
    *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{
      font-family: -apple-system, BlinkMacSystemFont, \'Segoe UI\', Roboto, sans-serif;
      background: #f8fafc;
      color: #1e293b;
      line-height: 1.6;
    }}
    .header {{
      background: #1e293b;
      color: #f8fafc;
      padding: 2rem 2.5rem;
    }}
    .header h1 {{
      font-size: 1.6rem;
      font-weight: 700;
      font-style: italic;
    }}
    .header .gca {{
      font-size: 0.95rem;
      color: #94a3b8;
      margin-top: 0.25rem;
      font-family: monospace;
    }}
    .badge {{
      display: inline-block;
      margin-top: 0.75rem;
      padding: 0.3rem 0.9rem;
      border-radius: 999px;
      font-size: 0.85rem;
      font-weight: 600;
      {badge_style};
    }}
    .container {{
      max-width: 960px;
      margin: 2rem auto;
      padding: 0 1.5rem;
    }}
    .cards {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
      gap: 1rem;
      margin-bottom: 2rem;
    }}
    .card {{
      background: #fff;
      border: 1px solid #e2e8f0;
      border-radius: 0.75rem;
      padding: 1.1rem 1.25rem;
      box-shadow: 0 1px 3px rgba(0,0,0,0.06);
    }}
    .card-label {{
      font-size: 0.75rem;
      font-weight: 600;
      text-transform: uppercase;
      letter-spacing: 0.05em;
      color: #64748b;
      margin-bottom: 0.4rem;
    }}
    .card-value {{
      font-size: 1.4rem;
      font-weight: 700;
      color: #0f172a;
    }}
    .card-sub {{
      font-size: 0.78rem;
      color: #94a3b8;
      margin-top: 0.2rem;
    }}
    .busco-bar-wrap {{
      background: #e2e8f0;
      border-radius: 999px;
      height: 10px;
      margin-top: 0.5rem;
      overflow: hidden;
    }}
    .busco-bar-fill {{
      height: 100%;
      border-radius: 999px;
    }}
    .protein-busco-bar-fill {{
      background: {protein_bar_color};
      width: {protein_bar_pct};
    }}
    .assembly-busco-bar-fill {{
      background: {assembly_bar_color};
      width: {assembly_bar_pct};
    }}
    .section {{
      background: #fff;
      border: 1px solid #e2e8f0;
      border-radius: 0.75rem;
      padding: 1.5rem;
      margin-bottom: 1.5rem;
      box-shadow: 0 1px 3px rgba(0,0,0,0.06);
    }}
    .section h2 {{
      font-size: 1rem;
      font-weight: 700;
      color: #1e293b;
      margin-bottom: 1rem;
      padding-bottom: 0.5rem;
      border-bottom: 1px solid #e2e8f0;
    }}
    .chart-container {{
      max-width: 600px;
      margin: 0 auto;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 0.9rem;
    }}
    tr:nth-child(even) {{ background: #f8fafc; }}
    td {{ padding: 0.55rem 0.75rem; border-bottom: 1px solid #f1f5f9; }}
    .metric-name {{
      font-weight: 600;
      color: #475569;
      width: 40%;
    }}
    .metric-value {{
      color: #1e293b;
      font-family: monospace;
      font-size: 0.88rem;
    }}
    .ftp-link {{
      color: #2563eb;
      word-break: break-all;
    }}
    .footer {{
      text-align: center;
      font-size: 0.8rem;
      color: #94a3b8;
      padding: 2rem 0;
    }}
  </style>
</head>
<body>

<div class="header">
  <h1>{html.escape(report.scientific_name)}</h1>
  <div class="gca">{html.escape(report.gca)}</div>
  <span class="badge">&#9679; {html.escape(report.protein_busco_quality)} Quality</span>
</div>

<div class="container">

  <!-- Key metric cards -->
  <div class="cards">
    <div class="card">
      <div class="card-label">Protein BUSCO</div>
      <div class="card-value">{protein_pct_display}</div>
      <div class="busco-bar-wrap"><div class="busco-bar-fill protein-busco-bar-fill"></div></div>
      <div class="card-sub">{html.escape(_fmt(report.protein_busco_lineage))}</div>
    </div>
    <div class="card">
      <div class="card-label">Assembly BUSCO</div>
      <div class="card-value">{assembly_pct_display}</div>
      <div class="busco-bar-wrap"><div class="busco-bar-fill assembly-busco-bar-fill"></div></div>
      <div class="card-sub">{html.escape(_fmt(report.assembly_busco_lineage))}</div>
    </div>
    <div class="card">
      <div class="card-label">Coding Genes</div>
      <div class="card-value">{html.escape(_fmt(report.coding_genes))}</div>
      <div class="card-sub">protein-coding</div>
    </div>
    <div class="card">
      <div class="card-label">Status</div>
      <div class="card-value" style="font-size:1rem">{html.escape(_fmt(report.gb_status))}</div>
      <div class="card-sub">{html.escape(_fmt(report.annotation_method))}</div>
    </div>
    <div class="card">
      <div class="card-label">Clade</div>
      <div class="card-value" style="font-size:1rem">
        {html.escape(_fmt(report.internal_clade))}
      </div>
      <div class="card-sub">Taxon ID: {html.escape(_fmt(report.lowest_taxon_id))}</div>
    </div>
    <div class="card">
      <div class="card-label">Release Date</div>
      <div class="card-value" style="font-size:1rem">{html.escape(_fmt(report.release_date))}</div>
      <div class="card-sub">Latest: {html.escape(_fmt(report.latest_annotated))}</div>
    </div>
  </div>

  <!-- BUSCO chart -->
  <div class="section">
    <h2>BUSCO Composition</h2>
    <div class="chart-container">
      <canvas id="buscoChart" height="160"></canvas>
    </div>
  </div>

  <!-- Full metrics table -->
  <div class="section">
    <h2>All Metrics</h2>
    <table>
      <tbody>
{table_rows_html}
        <tr>
          <td class="metric-name">FTP Link</td>
          <td class="metric-value">{ftp_html}</td>
        </tr>
      </tbody>
    </table>
  </div>

  <div class="footer">
    Generated by Ensembl Genebuild Metadata &mdash; {generated}
  </div>

</div>

<script>
  {chart_script}
</script>
</body>
</html>
"""

    with open(html_path, "w", encoding="utf-8") as f:
        f.write(html_doc)

    logger.info("Written HTML report to %s", html_path)
    return html_path
