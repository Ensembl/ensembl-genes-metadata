"""
clade_analysis.py

Module 2: PCA and MAD-based outlier detection across genome clades.

Takes the clade metrics DataFrame from clade_loader.load_clade_metrics()
and for each clade runs PCA to reduce the feature space, then uses
Median Absolute Deviation (MAD) to flag genomes whose annotation quality
is unusual compared to others in the same clade.

Why MAD instead of z-scores: MAD uses the median rather than the mean,
so it is not thrown off by the very outliers we are trying to detect.
Standard z-scores use the mean which gets skewed by extreme values.

Why per-clade PCA instead of one global PCA: feature distributions
differ too much across clades (e.g. fungi have very different gene
counts and intron lengths compared to vertebrates), so a global PCA
would mix biological differences with quality differences.

Note on BUSCO values: busco_completeness etc. are stored as raw gene
counts in the DB, not percentages. We convert them to percentages
here by dividing by the sum of complete + fragmented + missing.
"""

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from sklearn.decomposition import (
    PCA,
)  # pylint: disable=import-error  # type: ignore[import-untyped]  # noqa: E501
from sklearn.preprocessing import (
    StandardScaler,
)  # pylint: disable=import-error  # type: ignore[import-untyped]  # noqa: E501

logger = logging.getLogger(__name__)

# Minimum number of genomes needed in a clade to run PCA.
# Clades smaller than this get MAD on raw metrics only, no PCA.
MIN_PCA_SIZE = 10

# MAD outlier threshold, standard choice is 3.5 (Iglewicz and Hoaglin).
# A genome is flagged if its modified z-score exceeds this.
MAD_THRESHOLD = 8.0

# Features used for PCA, in order. These must exist as columns in the
# DataFrame returned by clade_loader after converting BUSCO to percentages.
_BUSCO_FEATURES = [
    "busco_completeness_pct",
    "busco_duplicated_pct",
    "busco_fragmented_pct",
    "busco_missing_pct",
]

_STATS_FEATURES = [
    "genebuild.stats.coding_genes",
    "genebuild.stats.total_transcripts",
    "genebuild.stats.transcripts_per_gene",
    "genebuild.stats.average_cds_length",
    "genebuild.stats.average_coding_intron_length",
    "genebuild.stats.single_exon_coding_genes",
    "genebuild.stats.overlapping_coding_genes",
]

ALL_FEATURES = _BUSCO_FEATURES + _STATS_FEATURES


@dataclass
class OutlierResult:  # pylint: disable=too-many-instance-attributes
    """Stores the outlier detection result for one genome in one clade."""

    gca: str
    scientific_name: str
    clade: str
    annotation_method: str
    is_outlier: bool
    mad_score: float
    outlier_features: List[str]
    pc1: Optional[float] = field(default=None)
    pc2: Optional[float] = field(default=None)
    busco_completeness_pct: Optional[float] = field(default=None)
    coding_genes: Optional[float] = field(default=None)
    clade_size: Optional[int] = field(default=None)


def _busco_to_pct(df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert raw BUSCO gene counts to percentages.

    The DB stores busco_completeness etc. as raw counts like 13098,
    not as percentages. We derive percentages by dividing each component
    by the total gene count (complete + fragmented + missing).
    """
    df = df.copy()

    for col in [
        "genebuild.busco_completeness",
        "genebuild.busco_duplicated",
        "genebuild.busco_fragmented",
        "genebuild.busco_missing",
        "genebuild.busco_single_copy",
    ]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    total = (
        df["genebuild.busco_completeness"]
        + df["genebuild.busco_fragmented"]
        + df["genebuild.busco_missing"]
    )
    total = total.replace(0, np.nan)

    df["busco_completeness_pct"] = df["genebuild.busco_completeness"] / total * 100
    df["busco_duplicated_pct"] = df["genebuild.busco_duplicated"] / total * 100
    df["busco_fragmented_pct"] = df["genebuild.busco_fragmented"] / total * 100
    df["busco_missing_pct"] = df["genebuild.busco_missing"] / total * 100

    return df


def _prepare_features(df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert all feature columns to numeric and drop rows where all
    stats features are null (genomes with no AGAT metrics).

    Rows missing only some features get those columns filled with
    the column median so PCA can still run on them.
    """
    df = df.copy()
    for col in _STATS_FEATURES:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    available = [f for f in ALL_FEATURES if f in df.columns]
    for col in available:
        median = df[col].median()
        df[col] = df[col].fillna(median)

    return df


def _mad_scores(values: np.ndarray) -> np.ndarray:
    """
    Compute modified MAD z-scores for a 1D array.

    Formula from Iglewicz and Hoaglin (1993):
    modified_z = 0.6745 * (x - median) / MAD
    where MAD = median(|x - median(x)|)

    Returns absolute values so we flag both high and low outliers.
    """
    median = np.median(values)
    mad = np.median(np.abs(values - median))
    if mad == 0:
        return np.zeros_like(values, dtype=float)
    return np.abs(0.6745 * (values - median) / mad)  # type: ignore[no-any-return]


def _run_pca(feature_matrix: np.ndarray, n_components: int = 2) -> np.ndarray:
    """
    Standardise and run PCA on a feature matrix.

    Returns the PCA-transformed coordinates (n_samples x n_components).
    Standardisation is needed because features have very different scales,
    e.g. busco_completeness_pct is 0-100 while coding_genes can be 20000+.
    """
    scaler = StandardScaler()
    scaled = scaler.fit_transform(feature_matrix)
    pca = PCA(n_components=min(n_components, feature_matrix.shape[1]))
    return pca.fit_transform(scaled)  # type: ignore[no-any-return]


def _analyse_clade(  # pylint: disable=too-many-locals
    clade_df: pd.DataFrame, clade_name: str
) -> List[OutlierResult]:
    """
    Run PCA and MAD outlier detection for one clade.

    For clades with enough genomes, runs PCA first and then MAD on the
    PC scores. Also runs MAD on each raw feature independently so we
    can report which specific metrics are unusual for each flagged genome.

    Args:
        clade_df: Subset of the full metrics DataFrame for this clade.
        clade_name: Name of the clade, used only for logging.

    Returns:
        List of OutlierResult, one per genome in the clade.
    """
    available_features = [f for f in ALL_FEATURES if f in clade_df.columns]
    if not available_features:
        logger.warning("No features available for clade %s, skipping", clade_name)
        return []

    feature_matrix = clade_df[available_features].values.astype(float)
    n_genomes = len(clade_df)

    # Run PCA if we have enough genomes
    pc_coords = None
    if n_genomes >= MIN_PCA_SIZE and len(available_features) >= 2:
        try:
            pc_coords = _run_pca(feature_matrix)
        except Exception as exc:  # pylint: disable=broad-except
            logger.warning("PCA failed for clade %s: %s", clade_name, exc)

    # MAD scores on each feature independently to find which ones are unusual
    feature_mad_scores = np.zeros((n_genomes, len(available_features)))
    for i, _ in enumerate(available_features):
        col_values = feature_matrix[:, i]
        feature_mad_scores[:, i] = _mad_scores(col_values)

    # Overall MAD score per genome: max across all features
    overall_mad = feature_mad_scores.max(axis=1)

    results = []
    for idx in range(n_genomes):
        row = clade_df.iloc[idx]
        is_outlier = bool(overall_mad[idx] > MAD_THRESHOLD)

        # Which specific features are driving the outlier flag
        outlier_features = [
            available_features[j]
            for j in range(len(available_features))
            if feature_mad_scores[idx, j] > MAD_THRESHOLD
        ]

        pc1 = float(pc_coords[idx, 0]) if pc_coords is not None else None
        pc2 = (
            float(pc_coords[idx, 1])
            if pc_coords is not None and pc_coords.shape[1] > 1
            else None
        )

        results.append(
            OutlierResult(
                gca=str(row["gca"]),
                scientific_name=str(row["scientific_name"]),
                clade=clade_name,
                annotation_method=str(row.get("annotation_method", "")),
                is_outlier=is_outlier,
                mad_score=float(overall_mad[idx]),
                outlier_features=outlier_features,
                pc1=pc1,
                pc2=pc2,
                clade_size=n_genomes,
                busco_completeness_pct=(
                    float(row["busco_completeness_pct"])
                    if "busco_completeness_pct" in row
                    and pd.notna(row["busco_completeness_pct"])
                    else None
                ),
                coding_genes=(
                    float(row["genebuild.stats.coding_genes"])
                    if "genebuild.stats.coding_genes" in row
                    and pd.notna(row["genebuild.stats.coding_genes"])
                    else None
                ),
            )
        )

    n_outliers = sum(1 for r in results if r.is_outlier)
    logger.info(
        "Clade %s: %d genomes, %d outliers flagged", clade_name, n_genomes, n_outliers
    )
    return results


def run_clade_analysis(
    clade_metrics_df: pd.DataFrame,
) -> Dict[str, List[OutlierResult]]:
    """
    Run PCA and MAD outlier detection across all clades.

    Args:
        clade_metrics_df: DataFrame from clade_loader.load_clade_metrics().

    Returns:
        Dict mapping clade name to list of OutlierResult for each genome
        in that clade. Clades with no usable features are skipped.
    """
    if clade_metrics_df.empty:
        logger.warning("Empty DataFrame passed to run_clade_analysis.")
        return {}

    df = _busco_to_pct(clade_metrics_df)
    df = _prepare_features(df)

    results: Dict[str, List[OutlierResult]] = {}
    for clade_name, clade_df in df.groupby("clade"):
        clade_results = _analyse_clade(clade_df.reset_index(drop=True), str(clade_name))
        if clade_results:
            results[str(clade_name)] = clade_results

    total_outliers = sum(
        1 for clade_results in results.values() for r in clade_results if r.is_outlier
    )
    logger.info(
        "Analysis complete: %d clades, %d total outliers flagged",
        len(results),
        total_outliers,
    )
    return results


def outlier_results_to_dataframe(
    results: Dict[str, List[OutlierResult]],
) -> pd.DataFrame:
    """
    Flatten the outlier results dict into a DataFrame for easy inspection
    and export to CSV.

    Args:
        results: Output of run_clade_analysis().

    Returns:
        DataFrame with one row per genome, sorted by clade then mad_score.
    """
    rows = []
    for clade_results in results.values():
        for r in clade_results:
            rows.append(
                {
                    "gca": r.gca,
                    "scientific_name": r.scientific_name,
                    "clade": r.clade,
                    "annotation_method": r.annotation_method,
                    "is_outlier": r.is_outlier,
                    "mad_score": round(r.mad_score, 3),
                    "outlier_features": ", ".join(r.outlier_features),
                    "pc1": r.pc1,
                    "pc2": r.pc2,
                    "busco_completeness_pct": r.busco_completeness_pct,
                    "coding_genes": r.coding_genes,
                    "clade_size": r.clade_size,
                }
            )
    if not rows:
        return pd.DataFrame()
    return (
        pd.DataFrame(rows)
        .sort_values(["clade", "mad_score"], ascending=[True, False])
        .reset_index(drop=True)
    )
