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
from typing import Optional

import pandas as pd
import pymysql  # pylint: disable=import-error
import pymysql.cursors  # pylint: disable=import-error

logger = logging.getLogger(__name__)

_METRICS_OF_INTEREST = [
    "genebuild.busco",
    "genebuild.busco_dataset",
    "genebuild.busco_version",
    "genebuild.stats.coding_genes",
]

_ANNO_WIDE_QUERY = """
SELECT
    CONCAT(a.gca_chain, '.', a.gca_version)   AS gca,
    s.scientific_name,
    s.common_name,
    a.lowest_taxon_id,
    s.clade                                      AS internal_clade,
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
    MAX(CASE WHEN am.metrics_name = 'genebuild.busco'
             THEN am.metrics_value END)           AS protein_busco,
    MAX(CASE WHEN am.metrics_name = 'genebuild.busco_dataset'
             THEN am.metrics_value END)           AS protein_busco_lineage,
    MAX(CASE WHEN am.metrics_name = 'genebuild.busco_version'
             THEN am.metrics_value END)           AS protein_busco_version,
    MAX(CASE WHEN am.metrics_name = 'genebuild.stats.coding_genes'
             THEN am.metrics_value END)           AS coding_genes
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
    AND am.metrics_name IN (
        'genebuild.busco',
        'genebuild.busco_dataset',
        'genebuild.busco_version',
        'genebuild.stats.coding_genes'
    )
GROUP BY
    a.assembly_id,
    gs.genebuild_status_id
"""


def _load_config(config_path: Optional[str] = None) -> dict:
    """Load database credentials from db_config.dev.json."""
    if config_path is None:
        config_path = str(
            Path(__file__).resolve().parents[6] / "db_config.dev.json"
        )
    with open(config_path, encoding="utf-8") as f:
        return json.load(f)


def get_connection(config_path: Optional[str] = None) -> pymysql.connections.Connection:
    """Return an open PyMySQL connection to gsoc_registry."""
    cfg = _load_config(config_path)
    conn = pymysql.connect(
        host=cfg["host"],
        port=int(cfg["port"]),
        user=cfg["user"],
        password=cfg.get("password", ""),
        database=cfg["database"],
        cursorclass=pymysql.cursors.DictCursor,
        connect_timeout=10,
    )
    logger.info("Connected to %s/%s", cfg["host"], cfg["database"])
    return conn


def load_anno_wide(
    config_path: Optional[str] = None,
    gca: Optional[str] = None,
) -> pd.DataFrame:
    """
    Query gsoc_registry and return a flat anno_wide DataFrame.

    Args:
        config_path: Path to db_config JSON. Defaults to db_config.dev.json
                     at the repo root.
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
