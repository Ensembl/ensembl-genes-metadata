"""
db_loader.py

Module 1: Database loader for gsoc_registry.

Connects to the Ensembl genebuild MySQL registry and returns an anno_wide
DataFrame compatible with extract_genome_report(). Uses PyMySQL with
credentials loaded from db_config.dev.json (gitignored).

The query pivots key-value rows from annotation_metrics into columns so
that the result matches the flat structure expected by GenomeReport.
"""

import json
import logging
from pathlib import Path
from typing import Dict, Optional

import pandas as pd
import pymysql  # pylint: disable=import-error
import pymysql.cursors  # pylint: disable=import-error

logger = logging.getLogger(__name__)

# Single source of truth: which annotation_metrics rows we pivot into
# columns, and what column alias each one becomes in anno_wide.
_METRICS_OF_INTEREST = [
    "genebuild.busco",
    "genebuild.busco_dataset",
    "genebuild.busco_version",
    "genebuild.stats.coding_genes",
]

_METRIC_TO_COLUMN: Dict[str, str] = {
    "genebuild.busco": "protein_busco",
    "genebuild.busco_dataset": "protein_busco_lineage",
    "genebuild.busco_version": "protein_busco_version",
    "genebuild.stats.coding_genes": "coding_genes",
}

# Assembly-level BUSCO lives in a SEPARATE table (assembly_metrics), not
# annotation_metrics. Confirmed via mentor review meeting transcript
# (Anna, 2026 review: "it's not in the annotation metrics, but the
# assembly metrics table -- we'll have to fetch from another table") and
# cross-checked against existing usage in annotations_service.py and
# assembly_service.py. Joined separately below via assembly_metrics.
_ASSEMBLY_METRICS_OF_INTEREST = [
    "assembly.busco",
    "assembly.busco_dataset",
    "assembly.busco_version",
]

_ASSEMBLY_METRIC_TO_COLUMN: Dict[str, str] = {
    "assembly.busco": "assembly_busco",
    "assembly.busco_dataset": "assembly_busco_lineage",
    "assembly.busco_version": "assembly_busco_version",
}

# Required keys for a valid db_config.dev.json. "password" is intentionally
# excluded since some local/dev configs use an empty/placeholder password.

# Metrics from the new_metrics table (AGAT-derived stats). These are stored
# separately from annotation_metrics and require their own join below.
_NEW_METRICS_OF_INTEREST = [
    "genebuild.stats.coding_genes",
    "genebuild.stats.total_transcripts",
    "genebuild.stats.transcripts_per_gene",
    "genebuild.stats.average_cds_length",
    "genebuild.stats.average_coding_intron_length",
    "genebuild.stats.single_exon_coding_genes",
    "genebuild.stats.longest_coding_gene_length",
    "genebuild.stats.average_coding_exon_length",
    "genebuild.stats.nc_non_coding_genes",
]
_NEW_METRIC_TO_COLUMN = {
    "genebuild.stats.coding_genes": "coding_genes",
    "genebuild.stats.total_transcripts": "total_transcripts",
    "genebuild.stats.transcripts_per_gene": "transcripts_per_gene",
    "genebuild.stats.average_cds_length": "average_cds_length",
    "genebuild.stats.average_coding_intron_length": "average_coding_intron_length",
    "genebuild.stats.single_exon_coding_genes": "single_exon_coding_genes",
    "genebuild.stats.longest_coding_gene_length": "longest_coding_gene_length",
    "genebuild.stats.average_coding_exon_length": "average_coding_exon_length",
    "genebuild.stats.nc_non_coding_genes": "nc_non_coding_genes",
}

_REQUIRED_CONFIG_KEYS = ("host", "port", "user", "database")


def _build_pivot_columns(
    metrics_list: list, column_map: Dict[str, str], table_alias: str
) -> str:
    """Build MAX(CASE WHEN ...) pivot columns for a given metrics table alias."""
    lines = []
    for metric_name in metrics_list:
        column_alias = column_map[metric_name]
        lines.append(
            f"    MAX(CASE WHEN {table_alias}.metrics_name = '{metric_name}'\n"
            f"             THEN {table_alias}.metrics_value END)           AS {column_alias}"
        )
    return ",\n".join(lines)


def _build_metrics_in_clause(metrics_list: list) -> str:
    """Build the SQL IN (...) clause from a list of metric names."""
    quoted = ", ".join(f"'{metric_name}'" for metric_name in metrics_list)
    return quoted


_ANNO_WIDE_QUERY = f"""
SELECT
    CONCAT(a.gca_chain, '.', a.gca_version)   AS gca,
    s.scientific_name,
    s.common_name,
    a.lowest_taxon_id,
    s.clade                                      AS internal_clade,
    gs.genebuild_status_id,
    gs.gb_status,
    gs.annotation_method,
    gs.genebuilder,
    gs.annotation_source,
    o.bioproject_id,
    gs.release_date,
    gs.date_status_update,
    gs.last_genebuild_update,
    gs.genebuild_version                         AS annotated_version,
    a.gca_version                                AS assembly_version,
{_build_pivot_columns(_METRICS_OF_INTEREST, _METRIC_TO_COLUMN, "am")},
{_build_pivot_columns(_ASSEMBLY_METRICS_OF_INTEREST, _ASSEMBLY_METRIC_TO_COLUMN, "asm")},
{_build_pivot_columns(_NEW_METRICS_OF_INTEREST, _NEW_METRIC_TO_COLUMN, "nm")}
FROM assembly a
JOIN species s
    ON s.lowest_taxon_id = a.lowest_taxon_id
JOIN genebuild_status gs
    ON gs.assembly_id = a.assembly_id
LEFT JOIN organism o
    ON o.assembly_id = a.assembly_id
LEFT JOIN annotation_metrics am
    ON am.assembly_id = a.assembly_id
    AND am.genebuild_status_id = gs.genebuild_status_id
    AND am.metrics_name IN ({_build_metrics_in_clause(_METRICS_OF_INTEREST)})
LEFT JOIN assembly_metrics asm
    ON asm.assembly_id = a.assembly_id
    AND asm.metrics_name IN ({_build_metrics_in_clause(_ASSEMBLY_METRICS_OF_INTEREST)})
LEFT JOIN new_metrics nm
    ON nm.assembly_id = a.assembly_id
    AND nm.genebuild_status_id = gs.genebuild_status_id
    AND nm.metrics_name IN ({_build_metrics_in_clause(_NEW_METRICS_OF_INTEREST)})
GROUP BY
    a.assembly_id,
    gs.genebuild_status_id
"""


def _load_config(config_path: Optional[str] = None) -> dict:
    """
    Load database credentials from db_config.dev.json.

    Validates that the loaded JSON is a dict and contains all required
    keys (host, port, user, database) so misconfiguration fails early
    and clearly, rather than surfacing as a confusing KeyError/TypeError
    deep inside pymysql.connect().
    """
    if config_path is None:
        config_path = str(
            Path(__file__).resolve().parents[5]
            / "backend"
            / "conf"
            / "db_config.dev.json"
        )

    with open(config_path, encoding="utf-8") as f:
        cfg = json.load(f)

    if not isinstance(cfg, dict):
        raise ValueError(
            f"Expected db_config to be a JSON object, got {type(cfg).__name__} "
            f"(path: {config_path})"
        )

    missing_keys = [key for key in _REQUIRED_CONFIG_KEYS if key not in cfg]
    if missing_keys:
        raise ValueError(
            f"db_config at {config_path} is missing required keys: {missing_keys}"
        )

    return cfg


def get_connection(config_path: Optional[str] = None) -> pymysql.connections.Connection:
    """Return an open PyMySQL connection to gsoc_registry."""
    cfg = _load_config(config_path)
    try:
        conn = pymysql.connect(
            host=cfg["host"],
            port=int(cfg["port"]),
            user=cfg["user"],
            password=cfg.get("password", ""),
            database=cfg["database"],
            cursorclass=pymysql.cursors.DictCursor,
            connect_timeout=10,
        )
    except pymysql.MySQLError as exc:
        logger.error(
            "Failed to connect to %s/%s: %s", cfg.get("host"), cfg.get("database"), exc
        )
        raise
    logger.info("Connected to %s/%s", cfg["host"], cfg["database"])
    return conn


def load_anno_wide(
    config_path: Optional[str] = None,
    gca: Optional[str] = None,
) -> pd.DataFrame:
    """
    Query gsoc_registry and return a flat anno_wide DataFrame.

    Args:
        config_path: Path to db_config JSON. Defaults to
                      metadata_app/backend/conf/db_config.dev.json.
        gca: Optional GCA accession to filter to a single genome.
             If None, returns all genomes.

    Returns:
        DataFrame with one row per genome annotation, columns matching
        the fields expected by extract_genome_report().
    """
    query = _ANNO_WIDE_QUERY
    params: tuple = ()
    if gca is not None:
        query += " HAVING gca = %s"
        params = (gca,)

    with get_connection(config_path) as conn:
        with conn.cursor() as cursor:
            cursor.execute(query, params)
            rows = cursor.fetchall()

    if not rows:
        logger.warning("No rows returned from gsoc_registry (gca filter: %s)", gca)
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    # Normalise date columns to strings so downstream code stays simple
    for col in ("release_date", "date_status_update", "last_genebuild_update"):
        if col in df.columns:
            df[col] = df[col].astype(str).replace("None", "")
    logger.info("Loaded %d rows from gsoc_registry", len(df))
    return df
