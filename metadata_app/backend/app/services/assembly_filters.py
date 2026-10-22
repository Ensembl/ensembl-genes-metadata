"""Shared assembly filtering rules used by assembly-facing services."""

import pandas as pd


CANDIDATE_MIN_CONTIG_N50 = 100_000
CANDIDATE_MIN_LOWEST_ILLUMINA_READS = 2
CANDIDATE_ASM_LEVELS = frozenset({"Complete genome", "Chromosome"})
CANDIDATE_UNTOUCHED_STATUS = "not_started"
CANDIDATE_REASSESSMENT_STATUSES = frozenset({"abandoned", "insufficient_data", "testing"})


def filter_by_genebuild_status(df, allowed_statuses):
    """Keep only GCAs whose status rows are all in ``allowed_statuses``."""
    status = df["gb_status"].fillna(CANDIDATE_UNTOUCHED_STATUS)
    blocked_gcas = set(df.loc[~status.isin(allowed_statuses), "gca"])
    return df.loc[~df["gca"].isin(blocked_gcas)].copy()


def filter_annotation_candidate_quality(df):
    """Keep assemblies matching candidate quality requirements, regardless of status."""
    required_columns = {
        "contig_n50",
        "short_read_paired_end_illumina_lowest",
        "asm_level",
        "is_current",
    }
    missing_columns = required_columns.difference(df.columns)
    if missing_columns:
        return df.iloc[0:0].copy()

    n50 = pd.to_numeric(df["contig_n50"], errors="coerce")
    lowest_illumina = pd.to_numeric(df["short_read_paired_end_illumina_lowest"], errors="coerce")

    quality_mask = (
        n50.gt(CANDIDATE_MIN_CONTIG_N50)
        & df["asm_level"].isin(CANDIDATE_ASM_LEVELS)
        & df["is_current"].eq("current")
        & lowest_illumina.ge(CANDIDATE_MIN_LOWEST_ILLUMINA_READS)
    )
    return df.loc[quality_mask].copy()


def filter_annotation_candidates(df, include_reassessment=False):
    """Apply the canonical annotation-candidate definition, including status."""
    if "gb_status" not in df.columns:
        return df.iloc[0:0].copy()

    quality_df = filter_annotation_candidate_quality(df)
    if quality_df.empty:
        return quality_df

    allowed_statuses = {CANDIDATE_UNTOUCHED_STATUS}
    if include_reassessment:
        allowed_statuses.update(CANDIDATE_REASSESSMENT_STATUSES)

    status = (
        df["gb_status"]
        .fillna(CANDIDATE_UNTOUCHED_STATUS)
        .replace({"not_annotated": CANDIDATE_UNTOUCHED_STATUS})
    )

    identity_column = "gca" if "gca" in df.columns else "assembly_id"
    if identity_column == "gca":
        status_mask = pd.Series(
            df.index.isin(filter_by_genebuild_status(df.assign(gb_status=status), allowed_statuses).index),
            index=df.index,
        )
    elif identity_column in df.columns:
        blocked_ids = set(df.loc[~status.isin(allowed_statuses), identity_column])
        status_mask = ~df[identity_column].isin(blocked_ids)
    else:
        status_mask = status.isin(allowed_statuses)

    return quality_df.loc[status_mask.reindex(quality_df.index, fill_value=False)].copy()
