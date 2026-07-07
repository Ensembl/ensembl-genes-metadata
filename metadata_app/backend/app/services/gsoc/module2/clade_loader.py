"""
clade_loader.py

Module 2: Clade-aware annotation metrics loader.

Loads annotation quality metrics from gsoc_registry grouped by biological
clade, using the taxonomy + clade_list tables to assign each genome to its
most specific matching clade. The resulting DataFrame is the input for
PCA and MAD-based outlier detection in clade_analysis.py.

Clade assignment strategy:
  - Each genome has a lowest_taxon_id in the species table.
  - The taxonomy table maps lowest_taxon_id to its full lineage
    (kingdom, phylum, class, order, family, genus, species), each as
    a taxon_class_id.
  - The clade_list table maps named clades (e.g. "aves", "mammalia")
    to their taxon_id.
  - A genome is assigned to a clade if any of its lineage taxon_class_ids
    matches a clade_list.taxon_id.
  - When multiple clades match (e.g. both "chordata" and "aves"), the
    most specific one (lowest count in clade_list, i.e. finest-grained)
    is used.

Note: species.clade is NULL for all rows in the current DB snapshot
(verified 2026-07-07). Clade assignment is therefore derived entirely
from the taxonomy + clade_list join.
"""

import logging
from typing import Optional

import pandas as pd
from metadata_app.backend.app.services.gsoc.module1.db_loader import (  # pylint: disable=import-error
    get_connection,
)

logger = logging.getLogger(__name__)

# Minimum number of genomes a clade must have (with at least one metric)
# before it is included in PCA/outlier analysis. Clades smaller than this
# do not have enough data points for meaningful statistics.
MIN_CLADE_SIZE = 10

# Metrics pulled from annotation_metrics (BUSCO) and new_metrics (AGAT stats)
# for use as PCA features. These were selected based on:
#   - Availability: present in >95% of genomes that have new_metrics rows
#   - Relevance: directly reflect annotation quality and completeness
#   - Non-redundancy: correlated metrics (e.g. average_cds_length vs
#     longest_cds_length) are not both included to avoid double-weighting
_PCA_METRICS = [
    "genebuild.busco",
    "genebuild.stats.coding_genes",
    "genebuild.stats.total_transcripts",
    "genebuild.stats.transcripts_per_gene",
    "genebuild.stats.average_cds_length",
    "genebuild.stats.average_coding_intron_length",
    "genebuild.stats.single_exon_coding_genes",
    "genebuild.stats.overlapping_coding_genes",
]

_CLADE_QUERY = """
SELECT
    CONCAT(a.gca_chain, '.', a.gca_version) AS gca,
    s.scientific_name,
    s.lowest_taxon_id,
    cl.clade,
    cl.taxon_rank,
    am.metrics_name,
    am.metrics_value
FROM assembly a
JOIN species s
    ON s.lowest_taxon_id = a.lowest_taxon_id
JOIN genebuild_status gs
    ON gs.assembly_id = a.assembly_id
    AND gs.gb_status = 'live'
JOIN taxonomy t
    ON t.lowest_taxon_id = s.lowest_taxon_id
JOIN clade_list cl
    ON cl.taxon_id = t.taxon_class_id
LEFT JOIN annotation_metrics am
    ON am.assembly_id = a.assembly_id
    AND am.genebuild_status_id = gs.genebuild_status_id
    AND am.metrics_name IN ({placeholders})
GROUP BY
    a.assembly_id,
    gs.genebuild_status_id,
    cl.clade_id,
    am.metrics_name
"""

_NEW_METRICS_QUERY = """
SELECT
    CONCAT(a.gca_chain, '.', a.gca_version) AS gca,
    nm.metrics_name,
    nm.metrics_value
FROM assembly a
JOIN genebuild_status gs
    ON gs.assembly_id = a.assembly_id
    AND gs.gb_status = 'live'
JOIN new_metrics nm
    ON nm.assembly_id = a.assembly_id
    AND nm.genebuild_status_id = gs.genebuild_status_id
    AND nm.metrics_name IN ({placeholders})
"""


def _build_placeholders(metrics: list) -> str:
    """Build a SQL IN clause placeholder string for a list of metric names."""
    return ", ".join(f"'{m}'" for m in metrics)


def _assign_finest_clade(df: pd.DataFrame) -> pd.DataFrame:
    """
    When a genome matches multiple clades, keep only the finest-grained one.

    Finest-grained = the clade with the fewest species in the DB (most
    specific taxonomic group). This is determined by counting how many
    distinct lowest_taxon_ids appear per clade in the loaded data, then
    keeping the clade with the smallest count for each GCA.

    Args:
        df: DataFrame with columns [gca, clade, ...], potentially with
            multiple clade rows per GCA.

    Returns:
        DataFrame with exactly one clade row per GCA.
    """
    clade_sizes = (
        df.groupby("clade")["gca"]
        .nunique()
        .reset_index()
        .rename(columns={"gca": "clade_size"})
    )
    df = df.merge(clade_sizes, on="clade", how="left")
    df = df.sort_values("clade_size")
    df = df.drop_duplicates(subset=["gca"], keep="first")
    df = df.drop(columns=["clade_size"])
    return df.reset_index(drop=True)


def load_clade_metrics(  # pylint: disable=too-many-locals
    config_path: Optional[str] = None,
) -> pd.DataFrame:
    """
    Load a wide-format DataFrame of annotation metrics with clade assignments.

    Each row represents one live-status genome with its clade assignment
    and one metric value. The result is pivoted by the caller into a
    genome x metric matrix for PCA.

    Args:
        config_path: Path to db_config JSON. Defaults to the standard
                     metadata_app config path.

    Returns:
        DataFrame with columns: gca, scientific_name, lowest_taxon_id,
        clade, taxon_rank, and one column per metric in _PCA_METRICS.
        Rows with no clade assignment are excluded.
        Clades with fewer than MIN_CLADE_SIZE genomes are excluded.
    """
    annotation_metrics = [
        m for m in _PCA_METRICS if not m.startswith("genebuild.stats")
    ]
    new_metrics_list = [m for m in _PCA_METRICS if m.startswith("genebuild.stats")]

    # Load BUSCO + clade assignment
    query = _CLADE_QUERY.format(placeholders=_build_placeholders(annotation_metrics))
    with get_connection(config_path) as conn:
        with conn.cursor() as cursor:
            cursor.execute(query)
            rows = cursor.fetchall()

    if not rows:
        logger.warning("No rows returned from clade metrics query.")
        return pd.DataFrame()

    base_df = pd.DataFrame(rows)
    logger.info("Loaded %d raw clade metric rows", len(base_df))

    # Pivot annotation_metrics into columns
    if "metrics_name" in base_df.columns and "metrics_value" in base_df.columns:
        identity_cols = [
            "gca",
            "scientific_name",
            "lowest_taxon_id",
            "clade",
            "taxon_rank",
        ]
        base_df = base_df.pivot_table(
            index=identity_cols,
            columns="metrics_name",
            values="metrics_value",
            aggfunc="first",
        ).reset_index()
        base_df.columns.name = None

    # Assign finest clade per genome
    base_df = _assign_finest_clade(base_df)

    # Load new_metrics (AGAT stats)
    if new_metrics_list:
        nm_query = _NEW_METRICS_QUERY.format(
            placeholders=_build_placeholders(new_metrics_list)
        )
        with get_connection(config_path) as conn:
            with conn.cursor() as cursor:
                cursor.execute(nm_query)
                nm_rows = cursor.fetchall()

        if nm_rows:
            nm_df = pd.DataFrame(nm_rows)
            nm_pivot = nm_df.pivot_table(
                index="gca",
                columns="metrics_name",
                values="metrics_value",
                aggfunc="first",
            ).reset_index()
            nm_pivot.columns.name = None
            base_df = base_df.merge(nm_pivot, on="gca", how="left")

    # Drop clades below minimum size
    clade_counts = base_df.groupby("clade")["gca"].nunique()
    valid_clades = clade_counts[clade_counts >= MIN_CLADE_SIZE].index
    before = len(base_df)
    base_df = base_df[base_df["clade"].isin(valid_clades)].reset_index(drop=True)
    dropped = before - len(base_df)
    if dropped:
        logger.info(
            "Dropped %d genomes from clades with fewer than %d members",
            dropped,
            MIN_CLADE_SIZE,
        )

    logger.info(
        "Final clade metrics DataFrame: %d genomes across %d clades",
        len(base_df),
        base_df["clade"].nunique() if not base_df.empty else 0,
    )
    return base_df
