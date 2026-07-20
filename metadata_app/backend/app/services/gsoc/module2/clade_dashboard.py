"""
clade_dashboard.py

Module 2: Clade summary dashboard generator.

Produces a single self-contained HTML file showing all 18 clades
in switchable tabs. Each tab shows:
  - Clade summary stats (genome count, outlier count, median BUSCO)
  - A sortable table of all genomes with their metrics and outlier status
  - Outliers highlighted in orange

Uses Chart.js for the BUSCO distribution chart per clade.
No external dependencies at render time except Chart.js from CDN.
"""

import html
import json
import logging
from datetime import date
from pathlib import Path


import pandas as pd

logger = logging.getLogger(__name__)


def _fmt(val: object) -> str:
    """Format a value for display, replacing None/NaN/empty with N/A."""
    if val is None:
        return "N/A"
    if isinstance(val, float):
        import math  # pylint: disable=import-outside-toplevel

        if math.isnan(val):
            return "N/A"
        return f"{val:.1f}"
    return str(val)


def _clade_stats(clade_df: pd.DataFrame) -> dict:
    """Compute summary statistics for one clade."""
    n_total = len(clade_df)
    n_outliers = int(clade_df["is_outlier"].sum())
    median_busco = clade_df["busco_completeness_pct"].median()
    median_coding = clade_df["coding_genes"].median()
    methods = clade_df["annotation_method"].value_counts().to_dict()
    return {
        "n_total": n_total,
        "n_outliers": n_outliers,
        "pct_outliers": round(n_outliers / n_total * 100, 1) if n_total else 0,
        "median_busco": (
            round(float(median_busco), 1) if not pd.isna(median_busco) else None
        ),
        "median_coding": (
            round(float(median_coding), 0) if not pd.isna(median_coding) else None
        ),
        "methods": methods,
    }


def _busco_histogram_data(clade_df: pd.DataFrame) -> str:
    """Return JSON array of BUSCO completeness values for Chart.js histogram."""
    values = clade_df["busco_completeness_pct"].dropna().tolist()
    bins = list(range(0, 101, 5))
    counts = [0] * len(bins)
    for v in values:
        idx = min(int(v // 5), len(bins) - 1)
        counts[idx] += 1
    labels = [f"{b}-{b+5}%" for b in bins]
    return json.dumps({"labels": labels, "counts": counts})


def _genome_rows_html(clade_df: pd.DataFrame) -> str:
    """Build HTML table rows for all genomes in a clade."""
    rows = []
    for _, row in clade_df.sort_values("mad_score", ascending=False).iterrows():
        is_outlier = bool(row["is_outlier"])
        bg = "background:#fff7ed" if is_outlier else ""
        outlier_badge = (
            "<span style='display:inline-block;padding:0.15rem 0.5rem;"
            "border-radius:999px;font-size:0.75rem;font-weight:600;"
            "color:#9a3412;background:#ffedd5'>Outlier</span>"
            if is_outlier
            else "<span style='color:#94a3b8;font-size:0.78rem'>-</span>"
        )
        rows.append(
            f"<tr style='border-bottom:1px solid #f1f5f9;{bg}'>"
            f"<td style='padding:0.4rem 0.6rem;font-family:monospace;font-size:0.78rem'>"
            f"{html.escape(str(row['gca']))}</td>"
            f"<td style='padding:0.4rem 0.6rem;font-size:0.82rem;font-style:italic'>"
            f"{html.escape(str(row['scientific_name']))}</td>"
            f"<td style='padding:0.4rem 0.6rem;font-size:0.82rem'>"
            f"{html.escape(str(row['annotation_method']))}</td>"
            f"<td style='padding:0.4rem 0.6rem;font-size:0.82rem;text-align:right'>"
            f"{_fmt(row['busco_completeness_pct'])}%</td>"
            f"<td style='padding:0.4rem 0.6rem;font-size:0.82rem;text-align:right'>"
            f"{_fmt(row['coding_genes'])}</td>"
            f"<td style='padding:0.4rem 0.6rem;font-size:0.82rem;text-align:right'>"
            f"{_fmt(row['mad_score'])}</td>"
            f"<td style='padding:0.4rem 0.6rem'>{outlier_badge}</td>"
            f"</tr>"
        )
    return "".join(rows)


def _clade_tab_content(_clade_name: str, clade_df: pd.DataFrame, idx: int) -> str:
    """Build the full HTML content for one clade tab."""
    stats = _clade_stats(clade_df)
    hist_data = _busco_histogram_data(clade_df)
    genome_rows = _genome_rows_html(clade_df)
    display = "block" if idx == 0 else "none"

    method_badges = " ".join(
        f"<span style='display:inline-block;padding:0.15rem 0.5rem;"
        f"border-radius:999px;font-size:0.75rem;background:#f1f5f9;color:#475569'>"
        f"{html.escape(m)}: {c}</span>"
        for m, c in stats["methods"].items()
    )

    busco_color = (
        "#16a34a"
        if (stats["median_busco"] or 0) >= 95
        else (
            "#d97706"
            if (stats["median_busco"] or 0) >= 85
            else "#ea580c" if (stats["median_busco"] or 0) >= 70 else "#dc2626"
        )
    )

    return f"""
<div id="tab-{idx}" class="tab-content" style="display:{display}">
  <div class="stats-grid">
    <div class="stat-card">
      <div class="stat-label">Total Genomes</div>
      <div class="stat-value">{stats["n_total"]}</div>
    </div>
    <div class="stat-card">
      <div class="stat-label">Outliers Flagged</div>
      <div class="stat-value" style="color:#9a3412">{stats["n_outliers"]}</div>
      <div class="stat-sub">{stats["pct_outliers"]}% of clade</div>
    </div>
    <div class="stat-card">
      <div class="stat-label">Median Protein BUSCO</div>
      <div class="stat-value" style="color:{busco_color}">
        {_fmt(stats["median_busco"])}%
      </div>
    </div>
    <div class="stat-card">
      <div class="stat-label">Median Coding Genes</div>
      <div class="stat-value">{_fmt(stats["median_coding"])}</div>
    </div>
  </div>

  <div class="section">
    <div style="font-size:0.82rem;color:#64748b;margin-bottom:0.75rem">
      Annotation methods: {method_badges}
    </div>
    <div style="max-width:600px;margin:0 auto">
      <canvas id="hist-{idx}" height="120"></canvas>
    </div>
    <script>
      new Chart(document.getElementById("hist-{idx}"), {{
        type: "bar",
        data: {{
          labels: {json.loads(hist_data)["labels"]},
          datasets: [{{
            label: "Genomes",
            data: {json.loads(hist_data)["counts"]},
            backgroundColor: "{busco_color}",
            borderRadius: 2,
          }}]
        }},
        options: {{
          responsive: true,
          plugins: {{
            legend: {{ display: false }},
            title: {{
              display: true,
              text: "BUSCO Completeness Distribution",
              font: {{ size: 12 }}
            }}
          }},
          scales: {{
            x: {{ title: {{ display: true, text: "BUSCO Completeness (%)" }} }},
            y: {{ title: {{ display: true, text: "Number of genomes" }} }}
          }}
        }}
      }});
    </script>
  </div>

  <div class="section">
    <h3 style="font-size:0.9rem;font-weight:700;color:#1e293b;
               margin-bottom:0.75rem">All Genomes</h3>
    <div style="overflow-x:auto">
      <table style="width:100%;border-collapse:collapse;font-size:0.85rem">
        <thead>
          <tr style="background:#f8fafc;border-bottom:2px solid #e2e8f0">
            <th style="padding:0.5rem 0.6rem;text-align:left;font-size:0.78rem;
                       color:#475569;text-transform:uppercase">GCA</th>
            <th style="padding:0.5rem 0.6rem;text-align:left;font-size:0.78rem;
                       color:#475569;text-transform:uppercase">Species</th>
            <th style="padding:0.5rem 0.6rem;text-align:left;font-size:0.78rem;
                       color:#475569;text-transform:uppercase">Method</th>
            <th style="padding:0.5rem 0.6rem;text-align:right;font-size:0.78rem;
                       color:#475569;text-transform:uppercase">BUSCO%</th>
            <th style="padding:0.5rem 0.6rem;text-align:right;font-size:0.78rem;
                       color:#475569;text-transform:uppercase">Coding Genes</th>
            <th style="padding:0.5rem 0.6rem;text-align:right;font-size:0.78rem;
                       color:#475569;text-transform:uppercase">MAD Score</th>
            <th style="padding:0.5rem 0.6rem;text-align:left;font-size:0.78rem;
                       color:#475569;text-transform:uppercase">Status</th>
          </tr>
        </thead>
        <tbody>{genome_rows}</tbody>
      </table>
    </div>
  </div>
</div>"""


def render_clade_dashboard(
    outlier_df: pd.DataFrame,
    output_path: Path,
) -> Path:
    """
    Render a single HTML clade summary dashboard with switchable tabs.

    Args:
        outlier_df: DataFrame from outlier_results_to_dataframe().
        output_path: Full path to write the HTML file.

    Returns:
        Path to the generated HTML file.
    """
    clades = sorted(outlier_df["clade"].unique())
    generated = date.today().isoformat()

    tab_buttons = []
    tab_contents = []

    for idx, clade_name in enumerate(clades):
        clade_df = outlier_df[
            outlier_df["clade"] == clade_name
        ].copy()  # pylint: disable=line-too-long
        n_outliers = int(clade_df["is_outlier"].sum())
        active = "active" if idx == 0 else ""
        tab_buttons.append(
            f"<button class='tab-btn {active}' onclick='showTab({idx})' id='btn-{idx}'>"
            f"{html.escape(clade_name.title())}"
            f"<span class='tab-count'>{len(clade_df)}</span>"
            f"{('<span class=tab-outlier>' + str(n_outliers) + ' outliers</span>') if n_outliers else ''}"  # pylint: disable=line-too-long
            f"</button>"
        )
        tab_contents.append(_clade_tab_content(clade_name, clade_df, idx))

    total_genomes = len(outlier_df)
    total_outliers = int(outlier_df["is_outlier"].sum())

    html_doc = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8"/>
  <meta name="viewport" content="width=device-width, initial-scale=1.0"/>
  <title>Clade Annotation Quality Dashboard</title>
  <script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
  <style>
    *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
      background: #f8fafc; color: #1e293b; line-height: 1.5;
    }}
    .header {{
      background: #1e293b; color: #f8fafc; padding: 1.5rem 2rem;
    }}
    .header h1 {{ font-size: 1.4rem; font-weight: 700; }}
    .header .sub {{
      font-size: 0.88rem; color: #94a3b8; margin-top: 0.25rem;
    }}
    .summary-bar {{
      background: #fff; border-bottom: 1px solid #e2e8f0;
      padding: 0.75rem 2rem; display: flex; gap: 2rem; align-items: center;
      font-size: 0.85rem; color: #475569;
    }}
    .summary-bar strong {{ color: #0f172a; }}
    .tab-bar {{
      background: #fff; border-bottom: 1px solid #e2e8f0;
      padding: 0 1.5rem; display: flex; flex-wrap: wrap; gap: 0.25rem;
      overflow-x: auto;
    }}
    .tab-btn {{
      padding: 0.6rem 1rem; border: none; background: none;
      cursor: pointer; font-size: 0.82rem; color: #64748b;
      border-bottom: 2px solid transparent; white-space: nowrap;
      display: flex; align-items: center; gap: 0.4rem;
    }}
    .tab-btn:hover {{ color: #1e293b; }}
    .tab-btn.active {{
      color: #1e293b; font-weight: 600;
      border-bottom: 2px solid #1e293b;
    }}
    .tab-count {{
      background: #f1f5f9; color: #64748b; padding: 0.1rem 0.4rem;
      border-radius: 999px; font-size: 0.72rem;
    }}
    .tab-outlier {{
      background: #ffedd5; color: #9a3412; padding: 0.1rem 0.4rem;
      border-radius: 999px; font-size: 0.72rem;
    }}
    .container {{ max-width: 1200px; margin: 1.5rem auto; padding: 0 1.5rem; }}
    .stats-grid {{
      display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
      gap: 1rem; margin-bottom: 1.5rem;
    }}
    .stat-card {{
      background: #fff; border: 1px solid #e2e8f0; border-radius: 0.75rem;
      padding: 1rem 1.25rem; box-shadow: 0 1px 3px rgba(0,0,0,0.05);
    }}
    .stat-label {{
      font-size: 0.72rem; font-weight: 600; text-transform: uppercase;
      letter-spacing: 0.05em; color: #64748b; margin-bottom: 0.3rem;
    }}
    .stat-value {{ font-size: 1.5rem; font-weight: 700; color: #0f172a; }}
    .stat-sub {{ font-size: 0.75rem; color: #94a3b8; margin-top: 0.15rem; }}
    .section {{
      background: #fff; border: 1px solid #e2e8f0; border-radius: 0.75rem;
      padding: 1.25rem; margin-bottom: 1.25rem;
      box-shadow: 0 1px 3px rgba(0,0,0,0.05);
    }}
    .footer {{
      text-align: center; font-size: 0.78rem; color: #94a3b8; padding: 2rem 0;
    }}
  </style>
</head>
<body>

<div class="header">
  <h1>Clade Annotation Quality Dashboard</h1>
  <div class="sub">
    Ensembl Genebuild Metadata &mdash; Module 2 Comparative Analysis
  </div>
</div>

<div class="summary-bar">
  <span><strong>{total_genomes}</strong> live genomes across
        <strong>{len(clades)}</strong> clades</span>
  <span><strong style="color:#9a3412">{total_outliers}</strong> outliers flagged
        ({round(total_outliers/total_genomes*100, 1)}% overall)</span>
  <span style="margin-left:auto;color:#94a3b8">Generated {generated}</span>
</div>

<div class="tab-bar">
  {"".join(tab_buttons)}
</div>

<div class="container">
  {"".join(tab_contents)}
  <div class="footer">
    Generated by Ensembl Genebuild Metadata &mdash; {generated}
  </div>
</div>

<script>
  function showTab(idx) {{
    document.querySelectorAll(".tab-content").forEach(function(el) {{
      el.style.display = "none";
    }});
    document.querySelectorAll(".tab-btn").forEach(function(el) {{
      el.classList.remove("active");
    }});
    document.getElementById("tab-" + idx).style.display = "block";
    document.getElementById("btn-" + idx).classList.add("active");
  }}
</script>

</body>
</html>
"""

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html_doc)

    logger.info("Written clade dashboard to %s", output_path)
    return output_path
