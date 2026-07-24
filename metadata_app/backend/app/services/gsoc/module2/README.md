# Module 2: Clade-Level Comparative Analysis and Outlier Detection

This is part of my GSoC 2026 project at EMBL-EBI, built on top of the Ensembl Assembly/Annotation Tracking App. Module 2 takes the per-genome metrics from Module 1 and puts them in context: for each genome, it asks whether its annotation quality looks unusual compared to other genomes in the same taxonomic group (clade). It uses PCA and MAD-based outlier detection to flag genomes that stand out.

The outlier results feed directly back into the Module 1 HTML reports, so when you open a per-genome report you can see at a glance whether that genome is unusual for its clade and which specific metrics are driving the flag.

## What this does, in plain terms

Knowing that a genome has 85% BUSCO completeness does not tell you much on its own. 85% might be excellent for a fungus and concerning for a primate. Module 2 answers the question: compared to other genomes in the same clade, is this annotation quality typical or unusual?

It does this by:
1. Loading annotation quality metrics for all live genomes from the registry, grouped by clade
2. Running PCA per clade to reduce the multi-metric feature space into a lower-dimensional representation
3. Using MAD (Median Absolute Deviation) to flag genomes whose metrics deviate significantly from the clade median
4. Returning outlier flags and scores that can be attached to individual genome reports

## Folder structure
## What each file actually does

**clade_loader.py**

This file handles everything to do with loading the data. It connects to the registry, pulls annotation metrics for all live genomes, and assigns each genome to its correct clade.

Clade assignment is the tricky part. The `species.clade` column in the database is NULL for all rows in the current snapshot, so we cannot use it directly. Instead, clade_loader uses the existing `taxonomy_service.assign_clade_and_species()` function together with `clade_settings.json` (the canonical clade definition file used by the genebuild pipelines) to walk each genome's taxonomic lineage and find its most specific matching clade.

The result is a wide DataFrame with one row per genome and one column per metric, ready for PCA.

A few things worth knowing:
- Only live-status genomes are included, for a clean consistent baseline
- Genomes with no matching clade are labelled Unassigned and excluded
- Clades with fewer than 10 genomes are excluded (not enough data for meaningful statistics)
- The legacy `clade_list` DB table is deliberately not used here. Anna confirmed it is legacy and will be deleted. `clade_settings.json` is the canonical source.

**clade_analysis.py**

This file takes the DataFrame from clade_loader and runs the actual analysis.

For each clade it:
1. Converts raw BUSCO gene counts to percentages (the DB stores them as raw counts like 13000, not as percentages)
2. Standardises all features using StandardScaler so metrics with different scales contribute equally to PCA
3. Runs PCA to get a lower-dimensional representation of each genome's quality profile
4. Computes MAD scores for each feature independently across the clade
5. Flags a genome as an outlier if its maximum MAD score exceeds the threshold (currently 8.0)
6. Returns an OutlierResult per genome with the clade name, clade size, MAD score, PC coordinates, and which specific features are driving the flag

Why MAD instead of z-scores: MAD uses the median rather than the mean, so it is not thrown off by the extreme values we are trying to detect. Standard z-scores can mask real outliers if the distribution is skewed.

Why per-clade PCA instead of one global PCA: feature distributions differ too much across clades. Fungi have very different gene counts and intron lengths compared to vertebrates, so a global PCA would conflate biological differences with annotation quality differences.

## PCA features used

The analysis uses 11 features per genome:

From annotation_metrics (BUSCO sub-metrics, converted to percentages):
- busco_completeness_pct
- busco_duplicated_pct
- busco_fragmented_pct
- busco_missing_pct

From new_metrics (AGAT-derived stats):
- genebuild.stats.coding_genes
- genebuild.stats.total_transcripts
- genebuild.stats.transcripts_per_gene
- genebuild.stats.average_cds_length
- genebuild.stats.average_coding_intron_length
- genebuild.stats.single_exon_coding_genes
- genebuild.stats.overlapping_coding_genes

## Clades currently covered

The pipeline loads 3,736 live genomes across 18 clades: ascomycota, aves, basidiomycota, cnidaria, coleoptera, diptera, hemiptera, hymenoptera, lepidoptera, lepidosauria, mammalia, mollusca, plants, porifera, rodentia, sharks, teleostei, testudines.

## How Module 2 connects to Module 1

The outlier results integrate into the per-genome HTML reports via `enrich_with_outlier_data()` in `genome_report.py`. Call it after `extract_genome_report()` to attach the outlier flag, MAD score, clade name, clade size, and flagged features to the report. The HTML renderer then shows a Clade Outlier Analysis section automatically if outlier data is present.

Example:

```python
from metadata_app.backend.app.services.gsoc.module1.db_loader import load_anno_wide
from metadata_app.backend.app.services.gsoc.module1.genome_report import (
    extract_genome_report,
    enrich_with_outlier_data,
)
from metadata_app.backend.app.services.gsoc.module1.html_renderer import render_html
from metadata_app.backend.app.services.gsoc.module2.clade_loader import load_clade_metrics
from metadata_app.backend.app.services.gsoc.module2.clade_analysis import run_clade_analysis
from pathlib import Path

# Run clade analysis once across all clades
clade_df = load_clade_metrics()
outlier_results = run_clade_analysis(clade_df)

# Generate a report for a specific genome enriched with outlier data
gca = "GCA_000002315.5"
anno_wide = load_anno_wide(gca=gca)
report = extract_genome_report(gca, anno_wide)
report = enrich_with_outlier_data(report, outlier_results)

# Render the HTML report with the outlier section
render_html(report, Path("outputs"))
```

## Setup

Same requirements as Module 1. See `module1/README.md` for the full setup instructions including the db_config.dev.json file and environment setup.

The database is only reachable from inside Codon or over the EBI VPN.

Additional dependency for Module 2:

```bash
pip install scikit-learn requests
```

## Running the tests

```bash
pytest metadata_app/backend/app/services/gsoc/module2/tests/ -v
```

37 tests currently pass across test_clade_loader.py and test_clade_analysis.py. All tests use synthetic fixture data and do not require a database connection.

## Code quality

All files are pylint 10/10, mypy clean, and black formatted. The only mypy suppressions are for sklearn and taxonomy_service imports which do not have type stubs available.

## Known limitations and things still to tune

**Outlier rate**: the current MAD threshold of 8.0 was chosen to reduce false positives, but some clades (particularly mammalia) still show high outlier rates because the clade clusters very tightly around 98-99% BUSCO completeness, making even minor deviations look extreme by MAD. This is under active discussion with mentors and the threshold or feature selection may be adjusted.

**Genomes missing from clade analysis**: some genomes in the registry do not have individual BUSCO sub-metric rows in annotation_metrics (only the composite string). These genomes are absent from the clade analysis pipeline because we need the numeric sub-metrics for PCA. This affects older annotations that predate the individual metric loading.

**Clade assignment for humans**: Homo sapiens is assigned to mammalia rather than primates because humans have a separate genebuild pipeline and are not in clade_settings.json as a distinct clade entry. This is expected and correct per Anna.

## MAD threshold — how it works and how to change it

The outlier detection uses a modified MAD (Median Absolute Deviation) score per feature. For each genome, the score is computed as:
The constant 0.6745 makes the MAD score equivalent to a standard z-score under a normal distribution. A genome is flagged as an outlier if its maximum MAD score across all features exceeds the threshold.

The current threshold is **8.0**, defined as `MAD_THRESHOLD` in `clade_analysis.py`:

```python
MAD_THRESHOLD = 8.0
```

To change it, update this constant. Lower values flag more genomes as outliers; higher values are more conservative. The value of 8.0 was chosen after observing that the default threshold of 3.5 produced false positive rates of 40-50% in some clades due to tight clustering. Some clades (particularly mammalia) still show high outlier rates at 8.0 because the clade clusters tightly around 98-99% BUSCO completeness, making even minor deviations appear extreme by MAD. Further tuning may be needed as more data is added to the registry.

## Notes on N/A fields in per-genome reports

Some fields in the per-genome HTML reports display N/A for certain genomes. These are genuine data gaps in the registry rather than code issues:

**FTP Link**: not all genomes have an FTP path stored in the registry. This field is populated from the existing anno_wide data and will show N/A if the registry does not have an FTP entry for that genome.

**Clade (card)**: the `species.clade` column is NULL for all rows in the current registry snapshot. Clade assignment is handled via `taxonomy_service` and `clade_settings.json` for Module 2 outlier detection, and the clade card in the per-genome report is populated from the outlier results when Module 2 is run. If a genome is not present in the clade analysis (e.g. it lacks individual BUSCO sub-metric rows), the clade card will remain N/A.

**Latest Annotated / Annotated Version**: these fields are NULL for many genomes in the registry and will show N/A where not populated.
