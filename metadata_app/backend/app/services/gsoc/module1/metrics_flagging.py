"""
metrics_flagging.py
Shared logic for classifying missing-metrics situations on a genebuild
status row.

A ``genebuild_status`` row can report a status (e.g. "completed", "live")
while having no corresponding rows in ``annotation_metrics`` and/or
``new_metrics``. Whether that is expected or worth investigating depends
entirely on the lifecycle stage (``gb_status``), confirmed directly with
the mentor team:

- Statuses where metrics are not expected yet (the build never reached a
  point where metrics would be computed) are not flagged at all.
- "completed" predates the introduction of automated metrics loading for
  a large fraction of rows, so a missing-metrics "completed" row is
  flagged as a legacy gap rather than an error.
- "check_busco" is sometimes assigned manually and can be messy; a
  missing-metrics row in this status is deliberately flagged so it forces
  reinvestigation.
- "archive" rows are old annotations that are no longer used and should
  be ignored completely - callers should exclude them from any
  comparison view entirely, not just suppress their flag.
- All other statuses (e.g. "live", "pre_released", "handed_over",
  "coming_soon") are expected to have metrics by that stage, so a missing
  row is flagged as needing investigation.

The same classification applies to both ``annotation_metrics`` and the
newer ``new_metrics`` table - call this function once per table with the
same ``gb_status`` and the relevant ``has_metrics`` value.

This module provides a single canonical classification location used by
both Module 1 (per-genome reports) and the multi-annotation comparison
view.
"""

from typing import Optional

from metadata_app.backend.app.services.gsoc.module1.logging_utils import (  # pylint: disable=import-error
    get_logger,
)

logger = get_logger(__name__)

# Canonical flag label vocabulary. This is the single source of truth for
# missing-metrics flag labels - renderers must source their labels/colours
# from these constants rather than redefining their own status checks or
# label strings.
FLAG_LEGACY_GAP = "Legacy Gap"
FLAG_INVESTIGATE = "Investigate"

# gb_status values where metrics are not expected at all (the build never
# reached a stage where metrics would be computed). No flag is shown.
_STATUSES_NO_METRICS_EXPECTED = frozenset(
    {"in_progress", "abandoned", "insufficient_data", "poor_genome_busco"}
)

# gb_status values where a missing-metrics row is a known, explained gap
# rather than something to investigate.
_STATUSES_LEGACY_GAP = frozenset({"completed"})

# gb_status values that should be excluded from missing-metrics
# evaluation, and from any comparison view, entirely.
_STATUSES_EXCLUDED = frozenset({"archive"})

# All other gb_status values (e.g. "live", "pre_released", "handed_over",
# "coming_soon", "check_busco") are expected to have metrics by that
# lifecycle stage, so a missing row is flagged as needing investigation.


def classify_missing_metrics(gb_status: str, has_metrics: bool) -> Optional[str]:
    """
    Classify whether a missing-metrics situation should be flagged.

    Args:
        gb_status: The genebuild_status.gb_status value for this row,
            e.g. "completed", "live", "check_busco".
        has_metrics: Whether at least one matching row exists in the
            metrics table being checked (annotation_metrics or
            new_metrics) for this genebuild_status row.

    Returns:
        FLAG_LEGACY_GAP if this is a known, explained gap (currently only
        "completed"); FLAG_INVESTIGATE if metrics are expected at this
        lifecycle stage but missing; None if no flag should be shown
        (metrics are present, not yet expected at this stage, or the
        status is "archive" and should be excluded entirely).

    Example:
        >>> classify_missing_metrics("completed", has_metrics=False)
        'Legacy Gap'
        >>> classify_missing_metrics("live", has_metrics=False)
        'Investigate'
        >>> classify_missing_metrics("live", has_metrics=True)
        >>> classify_missing_metrics("in_progress", has_metrics=False)
    """
    if not isinstance(gb_status, str) or not gb_status:
        logger.warning("Invalid gb_status received: %r", gb_status)
        return None

    if has_metrics:
        return None

    if gb_status in _STATUSES_EXCLUDED:
        return None

    if gb_status in _STATUSES_NO_METRICS_EXPECTED:
        return None

    if gb_status in _STATUSES_LEGACY_GAP:
        return FLAG_LEGACY_GAP

    return FLAG_INVESTIGATE


def is_excluded_status(gb_status: str) -> bool:
    """
    Check whether a genebuild_status row should be excluded entirely from
    comparison views, regardless of its metrics state.

    Currently only "archive" rows are excluded, per mentor confirmation
    that these are old annotations no longer in use.

    Args:
        gb_status: The genebuild_status.gb_status value for this row.

    Returns:
        True if the row should be excluded from comparison views entirely.

    Example:
        >>> is_excluded_status("archive")
        True
        >>> is_excluded_status("live")
        False
    """
    return gb_status in _STATUSES_EXCLUDED
