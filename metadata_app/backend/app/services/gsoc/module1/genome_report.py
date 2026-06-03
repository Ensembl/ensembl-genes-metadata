"""
genome_report.py

Module 1: Per-genome annotation metrics extraction.

Takes a single GCA accession, filters the anno_wide DataFrame produced by
the existing generate_tables() service, and returns a structured
GenomeReport dataclass ready for rendering into an HTML report.

This module does NOT query the database directly - it receives anno_wide
from the existing annotations_service.generate_tables() function, keeping
a clean separation between data retrieval and report generation.
"""

import logging
from dataclasses import dataclass, field
from typing import Optional
import pandas as pd

from metadata_app.backend.app.services.gsoc.module1.busco_utils import (  # pylint: disable=import-error
    parse_busco_string,
    busco_quality_label,
)


@dataclass
class GenomeReport:
    # pylint: disable=too-many-instance-attributes
    """  # pylint: disable=too-many-instance-attributes
    Structured representation of a single genome's annotation metrics.
    All fields map directly to columns in anno_wide from annotations_service.
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

    # Assembly BUSCO
    assembly_busco_raw: Optional[str]
    assembly_busco_complete: Optional[float]
    assembly_busco_lineage: Optional[str]
    assembly_busco_version: Optional[str]

    # Gene counts
    coding_genes: Optional[int]

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
    logging.info("Extracting genome report for GCA: %s", gca)

    genome_rows = anno_wide[anno_wide["gca"] == gca]

    if genome_rows.empty:
        logging.error("GCA %s not found in anno_wide", gca)
        raise ValueError(f"GCA '{gca}' not found in the provided dataset.")

    row = genome_rows.iloc[0]

    def safe_str(val) -> Optional[str]:
        if pd.isna(val) or val == "":
            return None
        return str(val)

    def safe_int(val) -> Optional[int]:
        try:
            if pd.isna(val) or val == "":
                return None
            return int(val)
        except (ValueError, TypeError):
            return None

    def safe_float(val) -> Optional[float]:
        try:
            if pd.isna(val) or val == "":
                return None
            return float(val)
        except (ValueError, TypeError):
            return None

    protein_busco_raw = safe_str(row.get("protein_busco"))
    protein_busco_parsed = parse_busco_string(protein_busco_raw or "")
    protein_busco_complete = protein_busco_parsed["complete"]

    assembly_busco_raw = safe_str(row.get("assembly_busco"))
    assembly_busco_parsed = parse_busco_string(assembly_busco_raw or "")
    assembly_busco_complete = assembly_busco_parsed["complete"]

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
        protein_busco_complete=(
            float(protein_busco_complete)
            if protein_busco_complete is not None
            else None
        ),
        protein_busco_lineage=safe_str(row.get("protein_busco_lineage")),
        protein_busco_version=safe_str(row.get("protein_busco_version")),
        protein_busco_quality=busco_quality_label(
            float(protein_busco_complete)
            if protein_busco_complete is not None
            else None
        ),
        assembly_busco_raw=assembly_busco_raw,
        assembly_busco_complete=(
            float(assembly_busco_complete)
            if assembly_busco_complete is not None
            else None
        ),
        assembly_busco_lineage=safe_str(row.get("assembly_busco_lineage")),
        assembly_busco_version=safe_str(row.get("assembly_busco_version")),
        coding_genes=safe_int(row.get("coding_genes")),
        latest_annotated=safe_str(row.get("latest_annotated")),
        annotated_version=safe_float(row.get("annotated_version")),
        assembly_version=safe_float(row.get("assembly_version")),
        ftp=safe_str(row.get("ftp")),
    )

    logging.info(
        "Successfully extracted report for %s (%s), BUSCO: %s%%, clade: %s",
        gca,
        report.scientific_name,
        report.protein_busco_complete,
        report.internal_clade,
    )
    return report
