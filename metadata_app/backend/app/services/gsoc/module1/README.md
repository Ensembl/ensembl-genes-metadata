# Module 1: Per-Genome Annotation Quality Reports

This is part of my GSoC 2026 project at EMBL-EBI, built on top of the Ensembl Assembly/Annotation Tracking App. Module 1 takes data from the `gsoc_registry` MySQL database and turns it into per-genome quality reports in four formats: HTML, CSV, TXT, and PNG plots.

## What this does, in plain terms

For a given genome assembly (identified by its GCA accession, e.g. `GCA_000001515.5`), this code pulls together its annotation metadata and BUSCO quality scores from the database, and produces a set of human-readable reports. The idea is that anyone, mentor or not, technical or not, should be able to open the HTML report or the PNG plots and immediately understand how good that genome's annotation is.

## Folder structure

```
metadata_app/backend/app/services/gsoc/module1/
    busco_utils.py        Parses raw BUSCO strings into numbers
    logging_utils.py      One shared logger setup used everywhere
    db_loader.py          Connects to gsoc_registry and pulls the data
    genome_report.py      Turns raw DB rows into a clean GenomeReport object
    report_renderer.py    Builds CSV, TXT, and PNG outputs
    html_renderer.py      Builds the HTML report
    tests/
        test_html_renderer.py
        test_report_renderer.py

tests/gsoc/module1/
    test_genome_report.py
```

## What each file actually does

**busco_utils.py**
BUSCO scores come out of the pipeline as a string like `C:97.7%[S:95.1%,D:2.6%],F:0.9%,M:1.4%,n:255`. This file is the one place in the whole codebase that knows how to read that string and pull the numbers out of it (completeness, single copy, duplicated, fragmented, missing, total gene count). It also decides what counts as "Excellent", "Good", "Moderate", or "Poor" quality, based on the completeness percentage. Every other file that needs to show a quality label or color gets it from here, so there's only one place to update if the thresholds ever change.

**logging_utils.py**
Just a small helper so every file logs messages the same way instead of each one setting up its own logger slightly differently.

**db_loader.py**
This is the only file that talks to the database. It connects to `gsoc_registry` using credentials from `db_config.dev.json` (this file is gitignored, see Setup below) and runs one query that joins together the assembly, species, genebuild status, and BUSCO metrics tables. It returns one row per genome with everything flattened into columns. Worth knowing: protein BUSCO comes from the `annotation_metrics` table, and assembly BUSCO comes from a separate table called `assembly_metrics`, joined independently. They look similar but they are not the same table.

**genome_report.py**
Takes the flat row from `db_loader.py` and turns it into a `GenomeReport`, which is a clean dataclass with one field per metric. This is also where we pick which annotation row to use if a genome has more than one (currently picks the most recently updated one, see Known Limitations below, this is a temporary rule pending mentor confirmation).

**report_renderer.py**
Takes a `GenomeReport` and writes out a CSV (all metrics in two columns), a TXT summary (readable plain text version), and two PNG plots (see below for what these mean).

**html_renderer.py**
Takes a `GenomeReport` and writes out a single self-contained HTML file. It uses Chart.js from a CDN for the BUSCO chart, so it needs an internet connection to render the chart properly, but the rest of the page works fine offline.

## What the PNG plots actually mean

Every genome report comes with two PNG images. Here's what to look at in each one.

**`<gca>_busco.png`** This is a side by side bar chart, one bar for Protein BUSCO and one for Assembly BUSCO. Each bar is split into four colored sections stacked on top of each other: green for single copy genes found correctly, dark green for duplicated genes, orange for fragmented genes, and red for genes that are missing entirely. A good annotation looks almost all green. If you see a lot of dark green (duplicated), that can mean the assembly has redundant or messy regions even if the overall completeness score looks high.

**`<gca>_quality_summary.png`** This is a quick at a glance card showing the genome's name, accession, both BUSCO percentages, the overall quality label, coding gene count, annotation method, status, clade, and whether it's the latest annotated version. The colored backgrounds match the same quality scale as the BUSCO chart, green is good, red is bad.

## Setup, for someone running this fresh

1. Clone the repo and check out this branch:
```bash
   git clone https://github.com/Shrinidhi-JN/ensembl-genes-metadata.git
   cd ensembl-genes-metadata
   git checkout shrinidhi/gsoc2026-annotation-metrics
```

2. Set up a Python 3.9 environment (Codon runs 3.9, so that's the version this code is written to support, even though local dev might use a newer version) and install dependencies:
```bash
   pip install -r requirements.txt --break-system-packages
```
   (or without `--break-system-packages` if you're using a virtual environment, which is recommended)

   Note: `requirements.txt` pins pandas to 2.3.3 specifically, not the newer 3.x line, because pandas 3.0+ requires Python 3.11 or higher and will not install on Codon's Python 3.9.

3. Create `metadata_app/backend/conf/db_config.dev.json` (this file is gitignored, you will not find it in the repo, you have to create it yourself):
```json
   {
     "host": "mysql-ens-genebuild-prod-1",
     "port": 4527,
     "user": "ensro",
     "database": "gsoc_registry",
     "password": ""
   }
```
   Note: this database is only reachable from inside Codon or over the EBI VPN. You will not be able to connect to it from a personal machine without one of those.

4. Run the tests to make sure everything is working:
```bash
   pytest tests/ metadata_app/backend/app/services/gsoc/module1/tests/ -v
```

5. To generate a report without needing a live database connection (useful for testing or demos), build a `GenomeReport` directly in Python and pass it to the renderers:
```python
   from metadata_app.backend.app.services.gsoc.module1.genome_report import GenomeReport
   from metadata_app.backend.app.services.gsoc.module1.report_renderer import render_report
   from metadata_app.backend.app.services.gsoc.module1.html_renderer import render_html

   report = GenomeReport(gca="GCA_000001515.5", ...)  # fill in all fields
   render_report(report, base_output_dir="outputs")
   render_html(report, output_dir)
```

6. To generate a report from the live database instead:
```python
   from metadata_app.backend.app.services.gsoc.module1.db_loader import load_anno_wide
   from metadata_app.backend.app.services.gsoc.module1.genome_report import extract_genome_report

   gca = "GCA_000001515.5"
   df = load_anno_wide(gca=gca)
   report = extract_genome_report(gca, df)
```
   Note: pass the whole DataFrame, not a single row. extract_genome_report() handles picking the right row internally (see Known Limitations below for how it currently does that when a GCA has more than one annotation).

## Code quality

Every file is run through pylint, mypy, and black before being committed. The one exception is a couple of pre-existing mypy warnings in `db_loader.py` about missing type stubs for the `pymysql` package itself, not anything in our code (these existed before this project started and are a repo wide thing, not specific to module1).

## Known limitations and things still being worked on

- **Multiple annotations per GCA**: right now if a genome has more than one annotation attempt in the database, we just pick the most recently updated one and log a warning. This is a stopgap, not the final design. The better options being discussed are either showing all annotations for a genome broken out by method, or generating a separate report per annotation. Not decided yet.
- **BUSCO sub-metrics**: the database actually stores duplication, fragmentation, and other BUSCO sub-scores as their own clean columns (not just buried inside the composite string we currently parse). We have not pivoted these in yet, but doing so would let us flag high duplication directly instead of relying on string parsing.
- **AGAT / new metrics**: gene/transcript/exon counts and length distributions from AGAT are not yet pulled into the reports. This needs a schema investigation session.

## Testing

104 tests currently pass across `tests/gsoc/module1/` and `metadata_app/backend/app/services/gsoc/module1/tests/`. Run them with:
```bash
pytest tests/ metadata_app/backend/app/services/gsoc/module1/tests/ -v
```
