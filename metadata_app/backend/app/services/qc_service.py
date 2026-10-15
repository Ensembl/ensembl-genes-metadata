"""Services for annotation quality-control data."""

import math
import re
import json
from pathlib import Path

import pandas as pd

from metadata_app.backend.app.core.database import get_db_connection
from metadata_app.backend.app.services.annotations_service import query_meta_registry

ANNOTATION_BUSCO_COLUMNS = {
    "protein_busco",
    "protein_busco_lineage",
    "protein_busco_version",
}

BUSCO_PATTERN = re.compile(
    r"C:(?P<complete>\d+(?:\.\d+)?)%.*?D:(?P<duplicated>\d+(?:\.\d+)?)%.*?"
    r"F:(?P<fragmented>\d+(?:\.\d+)?)%,M:(?P<missing>\d+(?:\.\d+)?)%"
)

STATUS_FILTERS = {
    "live": ["live"],
    "in progress": ["in_progress"],
    "handover": ["handed_over", "coming_soon"],
    "completed": ["completed", "pre_released"],
    "blocked": ["check_busco", "poor_genome_busco", "insufficient_data", "abandoned", "suppressed", "faulty", "large_sequence"],
}
ANNOTATION_METHOD_GROUPS = {
    "ensembl": ["projection_build", "full_genebuild", "anno"],
}
CLADE_SETTINGS_PATH = Path(__file__).resolve().parents[2] / "data" / "clade_settings.json"


def _status_label(raw_status):
    """Return the frontend status represented by a registry status."""
    if raw_status is None or pd.isna(raw_status):
        return raw_status
    for label, raw_statuses in STATUS_FILTERS.items():
        if str(raw_status) in raw_statuses:
            return label
    return raw_status


def _clade_taxon_id(clade: str):
    """Resolve a clade name from settings without making matching case-sensitive."""
    with CLADE_SETTINGS_PATH.open("r", encoding="utf-8") as handle:
        settings = json.load(handle)
    requested_clade = clade.casefold()
    for name, config in settings.items():
        if name.casefold() == requested_clade:
            return config.get("taxon_id")
    raise ValueError(f"Unknown QC clade: {clade}")


def _keep_qc_metrics(annotation_df: pd.DataFrame) -> pd.DataFrame:
    """Keep annotation BUSCO aliases and remove the other annotation metrics."""
    annotation_metric_columns = [
        column
        for column in annotation_df.columns
        if column.startswith("genebuild.") or column == "coding_genes"
    ]
    return annotation_df.drop(
        columns=[column for column in annotation_metric_columns if column not in ANNOTATION_BUSCO_COLUMNS]
    )


def _expand_annotation_methods(annotation_methods):
    """Expand frontend method groups into registry annotation methods."""
    if not annotation_methods:
        return None
    methods = [annotation_methods] if isinstance(annotation_methods, str) else annotation_methods
    expanded = []
    for method in methods:
        expanded.extend(ANNOTATION_METHOD_GROUPS.get(str(method).casefold(), [method]))
    return list(dict.fromkeys(expanded))


def _deduplicate_annotation_rows(annotation_df: pd.DataFrame) -> pd.DataFrame:
    """Keep one registry row per assembly version and annotation method."""
    identity = [column for column in ("gca", "annotation_method") if column in annotation_df.columns]
    if len(identity) < 2:
        return annotation_df
    return annotation_df.drop_duplicates(subset=identity, keep="first").reset_index(drop=True)


def query_new_metrics(annotation_df: pd.DataFrame) -> pd.DataFrame:
    """Return non-BUSCO new metrics for the annotations in ``annotation_df``."""
    if "genebuild_status_id" not in annotation_df.columns:
        raise ValueError("annotation_df must contain genebuild_status_id")

    genebuild_status_ids = sorted(annotation_df["genebuild_status_id"].dropna().unique().tolist())
    if not genebuild_status_ids:
        return pd.DataFrame(columns=["genebuild_status_id"])

    placeholders = ",".join(["%s"] * len(genebuild_status_ids))
    metrics_query = f"""
        SELECT genebuild_status_id, metrics_name, metrics_value
        FROM new_metrics
        WHERE genebuild_status_id IN ({placeholders})
          AND LOWER(metrics_name) NOT LIKE %s
    """

    with get_db_connection("meta") as conn:
        cursor = conn.cursor()
        cursor.execute(metrics_query, tuple(genebuild_status_ids) + ("%busco%",))
        metric_rows = cursor.fetchall()

    if not metric_rows:
        return pd.DataFrame(columns=["genebuild_status_id"])

    return (
        pd.DataFrame(metric_rows)
        .pivot_table(
            index="genebuild_status_id",
            columns="metrics_name",
            values="metrics_value",
            aggfunc="first",
        )
        .reset_index()
    )


def query_assembly_metrics(annotation_df: pd.DataFrame) -> pd.DataFrame:
    """Return assembly metrics needed for QC without changing annotations output."""
    if "assembly_id" not in annotation_df.columns:
        return pd.DataFrame(columns=["assembly_id"])

    assembly_ids = sorted(annotation_df["assembly_id"].dropna().unique().tolist())
    if not assembly_ids:
        return pd.DataFrame(columns=["assembly_id"])

    placeholders = ",".join(["%s"] * len(assembly_ids))
    metrics_query = f"""
        SELECT assembly_id, metrics_name, metrics_value
        FROM assembly_metrics
        WHERE assembly_id IN ({placeholders})
          AND metrics_name NOT IN (%s, %s, %s)
    """

    with get_db_connection("meta") as conn:
        cursor = conn.cursor()
        cursor.execute(
            metrics_query,
            tuple(assembly_ids)
            + ("assembly.busco", "assembly.busco_dataset", "assembly.busco_version"),
        )
        metric_rows = cursor.fetchall()

    if not metric_rows:
        return pd.DataFrame(columns=["assembly_id"])

    return (
        pd.DataFrame(metric_rows)
        .pivot_table(
            index="assembly_id",
            columns="metrics_name",
            values="metrics_value",
            aggfunc="first",
        )
        .reset_index()
    )


def query_qc_data(
    annotation_date,
    taxon_id,
    bioproject_id,
    group_name,
    gca=None,
    clade=None,
    status=None,
    annotation_method=None,
) -> pd.DataFrame:
    """Return annotation data plus QC-only metrics from ``new_metrics``."""
    if clade and not taxon_id:
        taxon_id = [_clade_taxon_id(clade)]

    raw_gcas = [gca] if isinstance(gca, str) else (gca or [])
    requested_gcas = [str(value).strip() for value in raw_gcas if str(value).strip()]
    taxonomy_cohort = bool(taxon_id) and not requested_gcas and not clade
    requested_annotation_methods = _expand_annotation_methods(annotation_method)
    if len(requested_gcas) == 1:
        requested_annotation_methods = _expand_annotation_methods("ensembl")
    annotation_df = query_meta_registry(
        annotation_date,
        taxon_id,
        bioproject_id,
        group_name,
        gca,
        include_internal_ids=True,
        gb_statuses=STATUS_FILTERS.get(status),
        annotation_methods=requested_annotation_methods,
    )
    annotation_df = _deduplicate_annotation_rows(annotation_df)

    # A single-GCA report still needs clade peers to make its outlier score
    # meaningful. Retrieve the requested record first, resolve its configured
    # clade, and use that clade as the scoring cohort. The API route filters
    # the prepared result back to the requested GCA before returning it.
    score_cohort = annotation_df
    if len(requested_gcas) == 1 and not annotation_df.empty:
        record_clade = annotation_df.iloc[0].get("internal_clade")
        if pd.notna(record_clade):
            try:
                score_cohort = query_meta_registry(
                    annotation_date,
                    [_clade_taxon_id(str(record_clade))],
                    bioproject_id,
                    group_name,
                    None,
                    include_internal_ids=True,
                    gb_statuses=["live"],
                    annotation_methods=requested_annotation_methods,
                )
                score_cohort = _deduplicate_annotation_rows(score_cohort)
                target_gca = str(annotation_df.iloc[0].get("gca", "")).casefold()
                cohort_gcas = {
                    str(value).casefold()
                    for value in score_cohort.get("gca", pd.Series(dtype=object)).dropna()
                }
                if target_gca and target_gca not in cohort_gcas:
                    score_cohort = pd.concat([score_cohort, annotation_df], ignore_index=True, sort=False)
            except (KeyError, ValueError):
                # Keep the requested record when its clade is not present in
                # the settings file; its score will correctly be 0 due to a
                # lack of a validated peer cohort.
                score_cohort = annotation_df

    annotation_df = score_cohort
    new_metrics_df = query_new_metrics(annotation_df)
    assembly_metrics_df = query_assembly_metrics(annotation_df)
    annotation_df = _keep_qc_metrics(annotation_df)
    new_metric_columns = [column for column in new_metrics_df.columns if column != "genebuild_status_id"]

    qc_df = annotation_df.merge(
        new_metrics_df,
        on="genebuild_status_id",
        how="left",
        validate="many_to_one",
    )
    qc_df = qc_df.merge(assembly_metrics_df, on="assembly_id", how="left", validate="many_to_one")
    qc_df.attrs["new_metric_columns"] = new_metric_columns
    qc_df.attrs["requested_gcas"] = requested_gcas
    qc_df.attrs["cohort_scope"] = "taxon" if taxonomy_cohort else "clade"
    return qc_df


def _first_value(row: pd.Series, names: tuple[str, ...]):
    """Find a metric by exact name or by its final dotted component."""
    for name in names:
        if name in row.index and pd.notna(row[name]):
            return row[name]
    for column in row.index:
        if column.rsplit(".", 1)[-1] in names and pd.notna(row[column]):
            return row[column]
    return None


def _number(value):
    if value is None or pd.isna(value):
        return None
    try:
        return float(str(value).replace(",", ""))
    except (TypeError, ValueError):
        return None


def _scaled_number(value, divisor=1.0, fraction_as_percent=False):
    number = _number(value)
    if number is None:
        return None
    if fraction_as_percent and 0 <= number <= 1:
        number *= 100
    return number / divisor


def _busco_values(value):
    if not value:
        return {}
    match = BUSCO_PATTERN.search(str(value))
    if not match:
        return {}
    return {key: float(number) for key, number in match.groupdict().items()}


def _sum_metrics(row: pd.Series, names: tuple[str, ...]):
    values = [_number(_first_value(row, (name,))) for name in names]
    values = [value for value in values if value is not None]
    return sum(values) if values else None


def _json_value(value):
    if value is None or pd.isna(value):
        return None
    if hasattr(value, "item"):
        value = value.item()
    if hasattr(value, "isoformat"):
        return value.isoformat()
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return None
    return value


def prepare_qc_records(qc_df: pd.DataFrame) -> list[dict]:
    """Shape registry rows into the stable record contract consumed by the QC page."""
    records = []
    for _, row in qc_df.iterrows():
        genome_busco = _busco_values(row.get("assembly_busco"))
        protein_busco = _busco_values(row.get("protein_busco"))
        new_metric_columns = qc_df.attrs.get("new_metric_columns", [])
        coding_genes = _number(_first_value(row, ("coding_genes",)))
        non_coding_genes = _number(_first_value(row, ("nc_non_coding_genes", "non_coding_genes")))
        pseudogenes = _number(_first_value(row, ("ps_pseudogenes", "pseudogenes")))

        records.append(
            {
                "gca": _json_value(row.get("gca")),
                "species": _json_value(row.get("scientific_name")) or "Unknown",
                "taxonomy": _json_value(row.get("internal_clade")) or "Unassigned",
                "taxonId": str(_json_value(row.get("lowest_taxon_id")) or ""),
                "status": _status_label(row.get("gb_status")),
                "associatedProject": _json_value(row.get("associated_project")),
                "groupName": _json_value(row.get("group_name")),
                "infraName": _json_value(row.get("infra_name")),
                "genebuilder": _json_value(row.get("genebuilder")),
                "annotationSource": _json_value(row.get("annotation_source")),
                "annotationMethod": _json_value(row.get("annotation_method")),
                "dateStarted": _json_value(row.get("date_started")),
                "releaseDate": _json_value(row.get("release_date")),
                "dateStatusUpdate": _json_value(row.get("date_status_update")),
                "lastGenebuildUpdate": _json_value(row.get("last_genebuild_update")),
                "assemblyBuscoLineage": _json_value(row.get("assembly_busco_lineage")),
                "assemblyBuscoVersion": _json_value(row.get("assembly_busco_version")),
                "contigN50Mb": _scaled_number(_first_value(row, ("contig_n50",)), 1_000_000),
                "scaffoldN50Mb": _scaled_number(_first_value(row, ("scaffold_n50",)), 1_000_000),
                "genomeSizeGb": _scaled_number(_first_value(row, ("genome_size", "genome_length")), 1_000_000_000),
                "gcPercent": _scaled_number(
                    _first_value(row, ("gc", "gc_content", "gc_percent")),
                    fraction_as_percent=True,
                ),
                "geneCount": sum(value for value in (coding_genes, non_coding_genes, pseudogenes) if value is not None)
                or None,
                "proteinCoding": coding_genes,
                "transcriptCount": _number(_first_value(row, ("total_transcripts", "transcript_count"))),
                "genomeBuscoComplete": genome_busco.get("complete"),
                "genomeBuscoDuplicated": genome_busco.get("duplicated"),
                "genomeBuscoFragmented": genome_busco.get("fragmented"),
                "genomeBuscoMissing": genome_busco.get("missing"),
                "proteinBuscoComplete": protein_busco.get("complete"),
                "proteinBuscoLineage": _json_value(row.get("protein_busco_lineage")),
                "proteinBuscoVersion": _json_value(row.get("protein_busco_version")),
                "buscoComplete": protein_busco.get("complete"),
                "buscoDuplicated": protein_busco.get("duplicated"),
                "buscoFragmented": protein_busco.get("fragmented"),
                "buscoMissing": protein_busco.get("missing"),
                "repeatPercent": _scaled_number(
                    _first_value(row, ("repeat_percent", "repeat_content", "repeats")),
                    fraction_as_percent=True,
                ),
                "annotationScore": protein_busco.get("complete"),
                "assemblyScore": genome_busco.get("complete"),
                "lncRnaCount": _number(_first_value(row, ("nc_long_non_coding_genes", "lnc_rna_count"))),
                "pseudogeneCount": pseudogenes,
                "smallRnaCount": _sum_metrics(
                    row,
                    (
                        "nc_mirna_genes",
                        "nc_rrna_genes",
                        "nc_scrna_genes",
                        "nc_snorna_genes",
                        "nc_snrna_genes",
                        "nc_trna_genes",
                        "nc_y_rna_genes",
                    ),
                ),
                "newMetrics": {
                    str(column): _json_value(row.get(column))
                    for column in new_metric_columns
                    if column in row.index and pd.notna(row.get(column))
                },
            }
        )

    for record in records:
        if qc_df.attrs.get("cohort_scope") == "taxon":
            peers = records
        else:
            peers = [
                peer
                for peer in records
                if peer["taxonomy"] == record["taxonomy"]
            ]
        features = ("proteinBuscoComplete", "genomeBuscoComplete", "geneCount", "transcriptCount")
        z_scores = []
        for feature in features:
            values = [peer[feature] for peer in peers if peer[feature] is not None]
            value = record[feature]
            if value is None or len(values) < 2:
                continue
            mean = sum(values) / len(values)
            variance = sum((item - mean) ** 2 for item in values) / len(values)
            standard_deviation = math.sqrt(variance)
            if standard_deviation:
                z_scores.append((value - mean) / standard_deviation)

        record["pcoa1"] = _zscore(record, peers, "proteinBuscoComplete")
        record["pcoa2"] = _zscore(record, peers, "geneCount")
        record["outlierScore"] = math.sqrt(sum(score * score for score in z_scores)) if z_scores else 0.0

    return records


def cohort_feature_means(records: list[dict], requested_gcas=None) -> dict[str, float | None]:
    """Return cohort means for the feature counts used by the gene-structure report."""
    feature_names = (
        "geneCount",
        "proteinCoding",
        "transcriptCount",
        "lncRnaCount",
        "pseudogeneCount",
        "smallRnaCount",
    )
    requested = {str(value).strip().casefold() for value in (requested_gcas or [])}
    peers = [record for record in records if str(record.get("gca", "")).casefold() not in requested]
    if not peers:
        peers = records
    return {
        feature: (
            sum(record[feature] for record in peers if record.get(feature) is not None)
            / sum(record.get(feature) is not None for record in peers)
            if any(record.get(feature) is not None for record in peers)
            else None
        )
        for feature in feature_names
    }


def _zscore(record: dict, peers: list[dict], field: str):
    value = record.get(field)
    values = [peer[field] for peer in peers if peer.get(field) is not None]
    if value is None or len(values) < 2:
        return 0.0
    mean = sum(values) / len(values)
    standard_deviation = math.sqrt(sum((item - mean) ** 2 for item in values) / len(values))
    return (value - mean) / standard_deviation if standard_deviation else 0.0
