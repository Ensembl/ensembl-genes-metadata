"""
genome_report.py

Module 1: Per-genome annotation metrics extraction.

Takes a single GCA accession, filters the anno_wide DataFrame produced by
the existing generate_tables() service, and returns a structured
GenomeReport dataclass ready for rendering into an HTML report.

This module does NOT query the database directly - it receives anno_wide
from the existing annotations_service.generate_tables() function, keeping
a clean separation between data retrieval and report generation.

TODO(AGAT/new-metrics): extend GenomeReport and the extraction below with
transcript/gene/exon counts and length-distribution summaries once the
new_metrics table mapping has been confirmed with the team (see Leanne's
2026-06-18 review feedback). A field-mapping table documenting source
table/metric name/type/missing-value handling should be added alongside.
"""

from dataclasses import dataclass, field
from typing import Any, Dict, Optional, Union

import pandas as pd

from metadata_app.backend.app.services.gsoc.module1.busco_utils import (  # pylint: disable=import-error
    parse_busco_string,
    busco_quality_label,
    busco_diff_label,
)
from metadata_app.backend.app.services.gsoc.module1.logging_utils import (  # pylint: disable=import-error
    get_logger,
)

logger = get_logger(__name__)


def safe_str(val: Any) -> Optional[str]:
    """Normalise a raw DataFrame value to a clean string, or None if missing."""
    if pd.isna(val) or val == "":
        return None
    return str(val)


def safe_int(val: Any) -> Optional[int]:
    """Normalise a raw DataFrame value to an int, or None if missing/invalid."""
    try:
        if pd.isna(val) or val == "":
            return None
        return int(val)
    except (ValueError, TypeError):
        return None


def safe_float(val: Any) -> Optional[float]:
    """Normalise a raw DataFrame value to a float, or None if missing/invalid."""
    try:
        if pd.isna(val) or val == "":
            return None
        return float(val)
    except (ValueError, TypeError):
        return None


def _busco_complete_as_float(
    parsed: Dict[str, Optional[Union[float, int, Dict[str, float]]]],
) -> Optional[float]:
    """
    Narrow parse_busco_string()["complete"] to Optional[float].

    parse_busco_string() returns a dict typed broadly enough to also hold
    the "extra" sub-dict, so callers that only want the numeric "complete"
    value need an explicit narrowing step rather than calling float()
    directly on a value mypy can't prove isn't a dict.
    """
    value = parsed["complete"]
    if isinstance(value, (float, int)):
        return float(value)
    return None


def _busco_extra_as_dict(
    parsed: Dict[str, Optional[Union[float, int, Dict[str, float]]]],
) -> Dict[str, float]:
    """Narrow parse_busco_string()["extra"] to a plain Dict[str, float]."""
    value = parsed["extra"]
    if isinstance(value, dict):
        return value
    return {}


@dataclass
class GenomeReport:  # pylint: disable=too-many-instance-attributes
    """
    Structured representation of a single genome's annotation metrics.
    All fields map directly to columns in anno_wide from annotations_service.

    Deliberately has many attributes (R0902 disabled above): this is a
    flat data-transfer object mirroring the columns of anno_wide, and
    splitting it into sub-objects would add indirection without benefit
    for a dataclass that is never extended with behaviour beyond field
    storage.
    """

    # Identity
    gca: str
    scientific_name: str
    common_name: Optional[str]
    lowest_taxon_id: Optional[int]
    internal_clade: Optional[str]

    # Annotation metadata
    annotation_method: Optional[str]
    genebuilder: Optional[str]
    gb_status: Optional[str]
    annotation_source: Optional[str]
    bioproject_id: Optional[str]
    associated_project: Optional[str]

    # Dates
    release_date: Optional[str]
    date_status_update: Optional[str]
    last_genebuild_update: Optional[str]

    # Protein BUSCO
    protein_busco_raw: Optional[str]
    protein_busco_complete: Optional[float]
    protein_busco_lineage: Optional[str]
    protein_busco_version: Optional[str]
    protein_busco_quality: str
    protein_busco_diff_label: str
    protein_busco_extra: Dict[str, float]

    # Assembly BUSCO
    assembly_busco_raw: Optional[str]
    assembly_busco_complete: Optional[float]
    assembly_busco_lineage: Optional[str]
    assembly_busco_version: Optional[str]
    assembly_busco_extra: Dict[str, float]

    # Gene counts and AGAT-derived stats from new_metrics
    coding_genes: Optional[int]
    total_transcripts: Optional[int]
    transcripts_per_gene: Optional[float]
    average_cds_length: Optional[float]
    average_coding_intron_length: Optional[float]
    single_exon_coding_genes: Optional[int]
    longest_coding_gene_length: Optional[int]
    average_coding_exon_length: Optional[float]
    nc_non_coding_genes: Optional[int]

    # Assembly version info
    latest_annotated: Optional[str]
    annotated_version: Optional[float]
    assembly_version: Optional[float]

    # FTP
    ftp: Optional[str]

    # Clade context (populated later by Module 2 if available)
    clade_median_busco: Optional[float] = field(default=None)
    clade_median_coding_genes: Optional[float] = field(default=None)
    clade_sample_size: Optional[int] = field(default=None)
    outlier_clade_size: Optional[int] = field(default=None)
    # Outlier detection results from Module 2 clade analysis
    is_outlier: Optional[bool] = field(default=None)
    outlier_mad_score: Optional[float] = field(default=None)
    outlier_features: Optional[str] = field(default=None)


def _select_current_annotation_row(genome_rows: pd.DataFrame, gca: str) -> pd.Series:
    """
    Pick the row representing the "current" annotation for a GCA when
    multiple genebuild_status rows exist for the same assembly.

    Selection rule (provisional - flagged for mentor confirmation, see
    code review 2026-06-18): the row with the most recent
    ``date_status_update`` wins. If that column is missing or all rows
    are unparseable, falls back to the first row and logs a warning so
    the ambiguity is visible rather than silently swallowed.

    Args:
        genome_rows: Rows from anno_wide already filtered to this GCA.
        gca: The GCA accession, used only for log messages.

    Returns:
        The single row (pandas Series) to use for this genome.
    """
    if len(genome_rows) == 1:
        return genome_rows.iloc[0]

    logger.warning(
        "Multiple annotation rows found for GCA %s (%d rows); "
        "selecting most recent by date_status_update.",
        gca,
        len(genome_rows),
    )

    if "date_status_update" not in genome_rows.columns:
        logger.warning(
            "date_status_update column unavailable for GCA %s; "
            "falling back to the first row.",
            gca,
        )
        return genome_rows.iloc[0]

    parsed_dates = pd.to_datetime(genome_rows["date_status_update"], errors="coerce")
    if parsed_dates.isna().all():
        logger.warning(
            "Could not parse date_status_update for any row of GCA %s; "
            "falling back to the first row.",
            gca,
        )
        return genome_rows.iloc[0]

    selected = genome_rows.loc[parsed_dates.idxmax()]
    assert isinstance(selected, pd.Series)
    return selected


def extract_genome_report(gca: str, anno_wide: pd.DataFrame) -> GenomeReport:
    """
    Extract a GenomeReport for a single GCA accession from anno_wide.

    Args:
        gca: GCA accession string e.g. "GCA_000001405.29"
        anno_wide: Full DataFrame returned by annotations_service.generate_tables()

    Returns:
        GenomeReport dataclass with all available metrics populated

    Raises:
        ValueError: If the GCA is not found in anno_wide
    """
    logger.info("Extracting genome report for GCA: %s", gca)

    genome_rows = anno_wide[anno_wide["gca"] == gca]

    if genome_rows.empty:
        logger.error("GCA %s not found in anno_wide", gca)
        raise ValueError(f"GCA '{gca}' not found in the provided dataset.")

    row = _select_current_annotation_row(genome_rows, gca)

    protein_busco_raw = safe_str(row.get("protein_busco"))
    protein_busco_parsed = parse_busco_string(protein_busco_raw or "")
    protein_busco_complete = _busco_complete_as_float(protein_busco_parsed)
    protein_busco_extra = _busco_extra_as_dict(protein_busco_parsed)

    assembly_busco_raw = safe_str(row.get("assembly_busco"))
    assembly_busco_parsed = parse_busco_string(assembly_busco_raw or "")
    assembly_busco_complete = _busco_complete_as_float(assembly_busco_parsed)
    assembly_busco_extra = _busco_extra_as_dict(assembly_busco_parsed)

    report = GenomeReport(
        gca=gca,
        scientific_name=safe_str(row.get("scientific_name")) or "Unknown",
        common_name=safe_str(row.get("common_name")),
        lowest_taxon_id=safe_int(row.get("lowest_taxon_id")),
        internal_clade=safe_str(row.get("internal_clade")),
        annotation_method=safe_str(row.get("annotation_method")),
        genebuilder=safe_str(row.get("genebuilder")),
        gb_status=safe_str(row.get("gb_status")),
        annotation_source=safe_str(row.get("annotation_source")),
        bioproject_id=safe_str(row.get("bioproject_id")),
        associated_project=safe_str(row.get("associated_project")),
        release_date=safe_str(row.get("release_date")),
        date_status_update=safe_str(row.get("date_status_update")),
        last_genebuild_update=safe_str(row.get("last_genebuild_update")),
        protein_busco_raw=protein_busco_raw,
        protein_busco_complete=protein_busco_complete,
        protein_busco_lineage=safe_str(row.get("protein_busco_lineage")),
        protein_busco_version=safe_str(row.get("protein_busco_version")),
        protein_busco_quality=busco_quality_label(protein_busco_complete),
        protein_busco_diff_label=busco_diff_label(
            protein_busco_complete, assembly_busco_complete
        ),
        protein_busco_extra=protein_busco_extra,
        assembly_busco_raw=assembly_busco_raw,
        assembly_busco_complete=assembly_busco_complete,
        assembly_busco_lineage=safe_str(row.get("assembly_busco_lineage")),
        assembly_busco_version=safe_str(row.get("assembly_busco_version")),
        assembly_busco_extra=assembly_busco_extra,
        coding_genes=safe_int(row.get("coding_genes")),
        total_transcripts=safe_int(row.get("total_transcripts")),
        transcripts_per_gene=safe_float(row.get("transcripts_per_gene")),
        average_cds_length=safe_float(row.get("average_cds_length")),
        average_coding_intron_length=safe_float(
            row.get("average_coding_intron_length")
        ),
        single_exon_coding_genes=safe_int(row.get("single_exon_coding_genes")),
        longest_coding_gene_length=safe_int(row.get("longest_coding_gene_length")),
        average_coding_exon_length=safe_float(row.get("average_coding_exon_length")),
        nc_non_coding_genes=safe_int(row.get("nc_non_coding_genes")),
        latest_annotated=safe_str(row.get("latest_annotated")),
        annotated_version=safe_float(row.get("annotated_version")),
        assembly_version=safe_float(row.get("assembly_version")),
        ftp=safe_str(row.get("ftp")),
    )

    logger.info(
        "Successfully extracted report for %s (%s), BUSCO: %s%%, clade: %s",
        gca,
        report.scientific_name,
        report.protein_busco_complete,
        report.internal_clade,
    )
    return report


def enrich_with_outlier_data(
    report: GenomeReport,
    outlier_results: "dict",
) -> GenomeReport:
    """
    Populate outlier fields on a GenomeReport from Module 2 clade analysis results.

    Args:
        report: A GenomeReport produced by extract_genome_report().
        outlier_results: Dict mapping clade name to list of OutlierResult,
                         as returned by clade_analysis.run_clade_analysis().

    Returns:
        The same GenomeReport with is_outlier, outlier_mad_score, and
        outlier_features populated if a matching result is found.
        Returns the report unchanged if no match is found.
    """
    if not outlier_results:
        return report

    # Search across all clades since internal_clade may be None
    # (species.clade is NULL in the current DB snapshot)
    match = None
    for clade_results in outlier_results.values():
        match = next((r for r in clade_results if r.gca == report.gca), None)
        if match:
            break

    if match is None:
        return report

    report.is_outlier = match.is_outlier
    report.outlier_mad_score = round(match.mad_score, 3)
    report.outlier_features = (
        ", ".join(match.outlier_features) if match.outlier_features else None
    )
    report.outlier_clade_size = match.clade_size
    if not report.internal_clade:
        report.internal_clade = match.clade

    logger.info(
        "Outlier data added for %s: is_outlier=%s, mad_score=%s, clade=%s, size=%s",
        report.gca,
        report.is_outlier,
        report.outlier_mad_score,
        report.internal_clade,
        report.outlier_clade_size,
    )
    return report
