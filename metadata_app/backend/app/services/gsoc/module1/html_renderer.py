"""
html_renderer.py

Module 1: HTML report generator.
Produces a self-contained HTML quality report from a GenomeReport.
No external dependencies at render time - Chart.js loaded from CDN.
"""

import logging
from datetime import date
from pathlib import Path

from metadata_app.backend.app.services.gsoc.module1.busco_utils import (  # pylint: disable=import-error
    parse_busco_string,
)
from metadata_app.backend.app.services.gsoc.module1.genome_report import (  # pylint: disable=import-error
    GenomeReport,
)

logger = logging.getLogger(__name__)


def _quality_badge_style(quality: str) -> str:
    """Return CSS colour pair for a quality label."""
    mapping = {
        "Excellent": ("#166534", "#dcfce7"),
        "High": ("#14532d", "#bbf7d0"),
        "Medium": ("#92400e", "#fef3c7"),
        "Low": ("#9a3412", "#ffedd5"),
        "Very Low": ("#7f1d1d", "#fee2e2"),
    }
    color, bg = mapping.get(quality, ("#374151", "#f3f4f6"))
    return f"color:{color};background:{bg}"


def _busco_bar_color(pct: float) -> str:
    """Return a hex fill colour for the BUSCO progress bar."""
    if pct >= 95:
        return "#16a34a"
    if pct >= 85:
        return "#d97706"
    if pct >= 70:
        return "#ea580c"
    return "#dc2626"


def _fmt(val: object) -> str:
    """Format a value for display, replacing None/empty with N/A."""
    if val is None or val == "":
        return "N/A"
    return str(val)


def _safe_float(v: object) -> float:
    """Safely cast a value to float, returning 0.0 on failure."""
    try:
        return float(v) if v is not None else 0.0  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return 0.0


def _busco_chart_js(report: GenomeReport) -> str:
    """Return an inline Chart.js script block for the BUSCO stacked bar."""
    p = parse_busco_string(report.protein_busco_raw or "")
    single = _safe_float(p.get("single_copy"))
    dup = _safe_float(p.get("duplicated"))
    frag = _safe_float(p.get("fragmented"))
    missing = _safe_float(p.get("missing"))

    ds_single = (
        f"{{ label: 'Single ({single:.1f}%)', data: [{single:.2f}],"
        f" backgroundColor: '#16a34a' }}"
    )
    ds_dup = (
        f"{{ label: 'Duplicated ({dup:.1f}%)', data: [{dup:.2f}],"
        f" backgroundColor: '#4ade80' }}"
    )
    ds_frag = (
        f"{{ label: 'Fragmented ({frag:.1f}%)', data: [{frag:.2f}],"
        f" backgroundColor: '#f59e0b' }}"
    )
    ds_missing = (
        f"{{ label: 'Missing ({missing:.1f}%)', data: [{missing:.2f}],"
        f" backgroundColor: '#ef4444' }}"
    )

    return f"""
    const ctx = document.getElementById('buscoChart').getContext('2d');
    new Chart(ctx, {{
        type: 'bar',
        data: {{
            labels: ['Protein BUSCO'],
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


def render_html(report: GenomeReport, output_dir: Path) -> Path:
    """
    Render a self-contained HTML quality report for a GenomeReport.

    Args:
        report: Populated GenomeReport dataclass.
        output_dir: Directory to write the HTML file into.

    Returns:
        Path to the generated HTML file.
    """
    html_path = output_dir / f"{report.gca.replace('.', '_')}_report.html"
    badge_style = _quality_badge_style(report.protein_busco_quality)
    busco_pct = report.protein_busco_complete or 0.0
    busco_bar_color = _busco_bar_color(busco_pct)
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
        ("Coding Genes", _fmt(report.coding_genes)),
        ("Latest Annotated", _fmt(report.latest_annotated)),
        ("Annotated Version", _fmt(report.annotated_version)),
        ("Assembly Version", _fmt(report.assembly_version)),
    ]

    table_rows_html = "\n".join(
        f"            <tr>\n"
        f'              <td class="metric-name">{name}</td>\n'
        f'              <td class="metric-value">{value}</td>\n'
        f"            </tr>"
        for name, value in metrics_rows
    )

    ftp_html = (
        f'<a href="{report.ftp}" target="_blank" class="ftp-link">' f"{report.ftp}</a>"
        if report.ftp
        else "N/A"
    )

    busco_bar_pct = f"{min(busco_pct, 100):.1f}%"

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8"/>
  <meta name="viewport" content="width=device-width, initial-scale=1.0"/>
  <title>Genome Report \u2014 {report.scientific_name} ({report.gca})</title>
  <script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
  <style>
    *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
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
      background: {busco_bar_color};
      width: {busco_bar_pct};
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
  <h1>{report.scientific_name}</h1>
  <div class="gca">{report.gca}</div>
  <span class="badge">&#9679; {report.protein_busco_quality} Quality</span>
</div>

<div class="container">

  <!-- Key metric cards -->
  <div class="cards">
    <div class="card">
      <div class="card-label">Protein BUSCO</div>
      <div class="card-value">{busco_pct:.1f}%</div>
      <div class="busco-bar-wrap"><div class="busco-bar-fill"></div></div>
      <div class="card-sub">{_fmt(report.protein_busco_lineage)}</div>
    </div>
    <div class="card">
      <div class="card-label">Coding Genes</div>
      <div class="card-value">{_fmt(report.coding_genes)}</div>
      <div class="card-sub">protein-coding</div>
    </div>
    <div class="card">
      <div class="card-label">Status</div>
      <div class="card-value" style="font-size:1rem">{_fmt(report.gb_status)}</div>
      <div class="card-sub">{_fmt(report.annotation_method)}</div>
    </div>
    <div class="card">
      <div class="card-label">Clade</div>
      <div class="card-value" style="font-size:1rem">{_fmt(report.internal_clade)}</div>
      <div class="card-sub">Taxon ID: {_fmt(report.lowest_taxon_id)}</div>
    </div>
    <div class="card">
      <div class="card-label">Release Date</div>
      <div class="card-value" style="font-size:1rem">{_fmt(report.release_date)}</div>
      <div class="card-sub">Latest: {_fmt(report.latest_annotated)}</div>
    </div>
  </div>

  <!-- BUSCO chart -->
  <div class="section">
    <h2>BUSCO Composition</h2>
    <div class="chart-container">
      <canvas id="buscoChart" height="120"></canvas>
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
        f.write(html)

    logger.info("Written HTML report to %s", html_path)
    return html_path
