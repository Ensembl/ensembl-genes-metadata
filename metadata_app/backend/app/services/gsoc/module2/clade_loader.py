"""
clade_loader.py

Module 2: Clade-aware annotation metrics loader.

Loads annotation quality metrics from gsoc_registry grouped by biological
clade, using the existing taxonomy_service.assign_clade_and_species()
function and clade_settings.json to assign each genome to its correct
clade. This replaces the previous clade_list DB table join, which Anna
confirmed is legacy and will be deleted.

Clade assignment strategy:
  - Load clade definitions from clade_settings.json (canonical source)
  - Load taxonomy hierarchy from the taxonomy DB table for all live genomes
  - For each genome, call assign_clade_and_species(lowest_taxon_id, ...)
    which walks species -> genus -> family -> order -> class -> phylum ->
    kingdom and returns the first matching clade from clade_settings.json
  - Genomes with no clade match are labelled "Unassigned" and excluded
  - Clades with fewer than MIN_CLADE_SIZE genomes are excluded from PCA

Note: Homo sapiens gets assigned to mammalia rather than primates because
humans have a separate pipeline and are not in clade_settings.json as a
distinct clade entry. This is expected and correct per Anna (2026-07-08).

Note: genebuild.busco is stored as a composite string in annotation_metrics
and cannot be used directly as a PCA feature. Individual numeric sub-metrics
(busco_completeness, busco_single_copy, etc.) are stored as separate rows
and are used instead.
"""

import logging
from typing import Dict, List, Optional

import pandas as pd

from metadata_app.backend.app.services.gsoc.module1.db_loader import (  # pylint: disable=import-error
    get_connection,
)
from metadata_app.backend.app.services.taxonomy_service import (  # pylint: disable=import-error
    assign_clade_and_species,
    load_clade_data,
)

logger = logging.getLogger(__name__)

# Minimum number of genomes a clade must have before it is included in
# PCA/outlier analysis. Clades smaller than this don't have enough data
# points for meaningful statistics.
MIN_CLADE_SIZE = 10

# Annotation metrics used as PCA features. These are individual numeric
# rows already stored separately in annotation_metrics and new_metrics.
# genebuild.busco is deliberately excluded because it is a composite
# string, not a number. The sub-metrics below are the numeric equivalents.
_ANNOTATION_PCA_METRICS = [
    "genebuild.busco_completeness",
    "genebuild.busco_single_copy",
    "genebuild.busco_duplicated",
    "genebuild.busco_fragmented",
    "genebuild.busco_missing",
]

_NEW_METRICS_PCA = [
    "genebuild.stats.coding_genes",
    "genebuild.stats.total_transcripts",
    "genebuild.stats.transcripts_per_gene",
    "genebuild.stats.average_cds_length",
    "genebuild.stats.average_coding_intron_length",
    "genebuild.stats.single_exon_coding_genes",
    "genebuild.stats.overlapping_coding_genes",
]

_TAXONOMY_QUERY = """
SELECT DISTINCT
    s.lowest_taxon_id,
    t.taxon_class_id,
    t.taxon_class
FROM species s
JOIN assembly a ON a.lowest_taxon_id = s.lowest_taxon_id
JOIN genebuild_status gs
    ON gs.assembly_id = a.assembly_id
    AND gs.gb_status = 'live'
JOIN taxonomy t ON t.lowest_taxon_id = s.lowest_taxon_id
"""

_ANNOTATION_METRICS_QUERY = """
SELECT
    CONCAT(a.gca_chain, '.', a.gca_version) AS gca,
    s.lowest_taxon_id,
    s.scientific_name,
    gs.annotation_method,
    am.metrics_name,
    am.metrics_value
FROM assembly a
JOIN species s ON s.lowest_taxon_id = a.lowest_taxon_id
JOIN genebuild_status gs
    ON gs.assembly_id = a.assembly_id
    AND gs.gb_status = 'live'
LEFT JOIN annotation_metrics am
    ON am.assembly_id = a.assembly_id
    AND am.genebuild_status_id = gs.genebuild_status_id
    AND am.metrics_name IN ({placeholders})
GROUP BY
    a.assembly_id,
    gs.genebuild_status_id,
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


def _build_placeholders(metrics: List[str]) -> str:
    """Build a SQL IN clause placeholder string for a list of metric names."""
    return ", ".join(f"'{m}'" for m in metrics)


def _load_taxonomy_dict(config_path: Optional[str]) -> Dict[str, list]:
    """
    Load the full taxonomy hierarchy from the DB for all live genomes.

    Returns a dict mapping str(lowest_taxon_id) to a list of
    {taxon_class, taxon_class_id} dicts, matching the format expected
    by taxonomy_service.assign_clade_and_species().
    """
    with get_connection(config_path) as conn:
        with conn.cursor() as cursor:
            cursor.execute(_TAXONOMY_QUERY)
            rows = cursor.fetchall()

    taxonomy_dict: Dict[str, list] = {}
    for row in rows:
        tid = str(row["lowest_taxon_id"])
        if tid not in taxonomy_dict:
            taxonomy_dict[tid] = []
        taxonomy_dict[tid].append(
            {
                "taxon_class": row["taxon_class"],
                "taxon_class_id": row["taxon_class_id"],
            }
        )

    logger.info("Loaded taxonomy hierarchy for %d taxa", len(taxonomy_dict))
    return taxonomy_dict


def _load_annotation_metrics(config_path: Optional[str]) -> pd.DataFrame:
    """
    Load annotation metrics (BUSCO sub-metrics) for all live genomes.

    Returns a wide DataFrame with one row per genome and one column
    per metric in _ANNOTATION_PCA_METRICS.
    """
    query = _ANNOTATION_METRICS_QUERY.format(
        placeholders=_build_placeholders(_ANNOTATION_PCA_METRICS)
    )
    with get_connection(config_path) as conn:
        with conn.cursor() as cursor:
            cursor.execute(query)
            rows = cursor.fetchall()

    if not rows:
        logger.warning("No annotation metrics rows returned.")
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    identity_cols = ["gca", "lowest_taxon_id", "scientific_name", "annotation_method"]
    pivot = df.pivot_table(
        index=identity_cols,
        columns="metrics_name",
        values="metrics_value",
        aggfunc="first",
    ).reset_index()
    pivot.columns.name = None
    return pivot


def _load_new_metrics(config_path: Optional[str]) -> pd.DataFrame:
    """
    Load AGAT new_metrics for all live genomes.

    Returns a wide DataFrame with one row per genome and one column
    per metric in _NEW_METRICS_PCA.
    """
    query = _NEW_METRICS_QUERY.format(
        placeholders=_build_placeholders(_NEW_METRICS_PCA)
    )
    with get_connection(config_path) as conn:
        with conn.cursor() as cursor:
            cursor.execute(query)
            rows = cursor.fetchall()

    if not rows:
        logger.warning("No new_metrics rows returned.")
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    pivot = df.pivot_table(
        index="gca",
        columns="metrics_name",
        values="metrics_value",
        aggfunc="first",
    ).reset_index()
    pivot.columns.name = None
    return pivot


def _assign_clades(
    base_df: pd.DataFrame,
    taxonomy_dict: Dict[str, list],
    clade_data: dict,
) -> pd.DataFrame:
    """
    Assign clade to each genome using taxonomy_service.

    Adds a clade and annotation_method column to base_df.
    Genomes with no clade match get "Unassigned".
    """
    clades = []
    for _, row in base_df.iterrows():
        taxon_id = row["lowest_taxon_id"]
        clade, _, _ = assign_clade_and_species(taxon_id, clade_data, taxonomy_dict)
        clades.append(clade)
    base_df = base_df.copy()
    base_df["clade"] = clades
    return base_df


def _drop_small_clades(base_df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove genomes belonging to clades with fewer than MIN_CLADE_SIZE members.

    Args:
        base_df: DataFrame with a clade column.

    Returns:
        Filtered DataFrame with small clades and Unassigned removed.
    """
    before = len(base_df)
    base_df = base_df[base_df["clade"] != "Unassigned"]
    clade_counts = base_df.groupby("clade")["gca"].nunique()
    valid_clades = clade_counts[clade_counts >= MIN_CLADE_SIZE].index
    result = base_df[base_df["clade"].isin(valid_clades)].reset_index(drop=True)
    dropped = before - len(result)
    if dropped:
        logger.info(
            "Dropped %d genomes (Unassigned or clades < %d members)",
            dropped,
            MIN_CLADE_SIZE,
        )
    return result


def load_clade_metrics(  # pylint: disable=too-many-locals
    config_path: Optional[str] = None,
) -> pd.DataFrame:
    """
    Load a wide-format DataFrame of annotation metrics with clade assignments.

    Uses taxonomy_service.assign_clade_and_species() and clade_settings.json
    (the canonical clade source) rather than the legacy clade_list DB table.

    Args:
        config_path: Path to db_config JSON. Defaults to the standard
                     metadata_app config path.

    Returns:
        DataFrame with columns: gca, scientific_name, lowest_taxon_id,
        annotation_method, clade, and one column per PCA metric.
        Unassigned genomes and clades below MIN_CLADE_SIZE are excluded.
    """
    logger.info("Loading clade settings from clade_settings.json")
    clade_data = load_clade_data()

    logger.info("Loading taxonomy hierarchy from DB")
    taxonomy_dict = _load_taxonomy_dict(config_path)

    logger.info("Loading annotation metrics from DB")
    base_df = _load_annotation_metrics(config_path)
    if base_df.empty:
        return pd.DataFrame()

    logger.info("Loading new_metrics from DB")
    nm_df = _load_new_metrics(config_path)
    if not nm_df.empty:
        base_df = base_df.merge(nm_df, on="gca", how="left")

    logger.info("Assigning clades via taxonomy_service")
    base_df = _assign_clades(base_df, taxonomy_dict, clade_data)

    base_df = _drop_small_clades(base_df)

    logger.info(
        "Final clade metrics DataFrame: %d genomes across %d clades",
        len(base_df),
        base_df["clade"].nunique() if not base_df.empty else 0,
    )
    return base_df
