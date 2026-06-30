"""
test_metrics_flagging.py

Unit tests for metrics_flagging.py.
All tests are pure-function checks against the gb_status flagging rules
confirmed with the mentor team - no database connection required.
"""

from metadata_app.backend.app.services.gsoc.module1.metrics_flagging import (
    FLAG_INVESTIGATE,
    FLAG_LEGACY_GAP,
    classify_missing_metrics,
    is_excluded_status,
)

# classify_missing_metrics - statuses where metrics are not expected at all


def test_in_progress_missing_not_flagged() -> None:
    """in_progress with no metrics is not flagged - build is not finished."""
    assert classify_missing_metrics("in_progress", has_metrics=False) is None


def test_abandoned_missing_not_flagged() -> None:
    """abandoned with no metrics is not flagged - work was stopped."""
    assert classify_missing_metrics("abandoned", has_metrics=False) is None


def test_insufficient_data_missing_not_flagged() -> None:
    """insufficient_data with no metrics is not flagged - never had enough
    data to attempt the build."""
    assert classify_missing_metrics("insufficient_data", has_metrics=False) is None


def test_poor_genome_busco_missing_not_flagged() -> None:
    """poor_genome_busco with no metrics is not flagged - failed QC before
    annotation metrics would be computed."""
    assert classify_missing_metrics("poor_genome_busco", has_metrics=False) is None


# classify_missing_metrics - completed (legacy gap)


def test_completed_missing_is_legacy_gap() -> None:
    """completed with no metrics is flagged as a known legacy gap, per
    mentor confirmation that completed annotations come from an old
    pipeline that predates automated metrics loading."""
    assert classify_missing_metrics("completed", has_metrics=False) == FLAG_LEGACY_GAP


def test_completed_with_metrics_not_flagged() -> None:
    """completed with metrics present is not flagged at all."""
    assert classify_missing_metrics("completed", has_metrics=True) is None


# classify_missing_metrics - check_busco (deliberately flagged)


def test_check_busco_missing_is_investigate() -> None:
    """check_busco with no metrics is flagged as needing investigation,
    per mentor confirmation that this status is manually assigned and can
    get messy - the flag is meant to force reinvestigation."""
    assert (
        classify_missing_metrics("check_busco", has_metrics=False) == FLAG_INVESTIGATE
    )


def test_check_busco_with_metrics_not_flagged() -> None:
    """check_busco with metrics present is not flagged at all."""
    assert classify_missing_metrics("check_busco", has_metrics=True) is None


# classify_missing_metrics - archive (excluded entirely)


def test_archive_missing_not_flagged() -> None:
    """archive with no metrics is not flagged - these rows should be
    ignored completely, per mentor confirmation that archived annotations
    are old and no longer used."""
    assert classify_missing_metrics("archive", has_metrics=False) is None


def test_archive_with_metrics_not_flagged() -> None:
    """archive with metrics present is also not flagged - archive status
    is excluded from this logic regardless of metrics presence."""
    assert classify_missing_metrics("archive", has_metrics=True) is None


# classify_missing_metrics - statuses where metrics are expected


def test_live_missing_is_investigate() -> None:
    """live with no metrics is flagged as needing investigation - live
    annotations are expected to have metrics by this stage."""
    assert classify_missing_metrics("live", has_metrics=False) == FLAG_INVESTIGATE


def test_live_with_metrics_not_flagged() -> None:
    """live with metrics present is not flagged at all."""
    assert classify_missing_metrics("live", has_metrics=True) is None


def test_pre_released_missing_is_investigate() -> None:
    """pre_released with no metrics is flagged as needing investigation."""
    assert (
        classify_missing_metrics("pre_released", has_metrics=False) == FLAG_INVESTIGATE
    )


def test_handed_over_missing_is_investigate() -> None:
    """handed_over with no metrics is flagged as needing investigation."""
    assert (
        classify_missing_metrics("handed_over", has_metrics=False) == FLAG_INVESTIGATE
    )


def test_coming_soon_missing_is_investigate() -> None:
    """coming_soon with no metrics is flagged as needing investigation."""
    assert (
        classify_missing_metrics("coming_soon", has_metrics=False) == FLAG_INVESTIGATE
    )


# classify_missing_metrics - invalid input handling


def test_none_status_returns_none() -> None:
    """A None gb_status does not crash and returns no flag."""
    result = classify_missing_metrics(None, has_metrics=False)  # type: ignore[arg-type]
    assert result is None


def test_empty_string_status_returns_none() -> None:
    """An empty-string gb_status does not crash and returns no flag."""
    assert classify_missing_metrics("", has_metrics=False) is None


def test_unrecognised_status_defaults_to_investigate() -> None:
    """A gb_status not in any known bucket defaults to investigate rather
    than silently passing, so unexpected/new enum values surface instead
    of being missed."""
    assert (
        classify_missing_metrics("some_future_status", has_metrics=False)
        == FLAG_INVESTIGATE
    )


# is_excluded_status


def test_is_excluded_status_archive_is_true() -> None:
    """archive is the only status excluded from comparison views entirely."""
    assert is_excluded_status("archive") is True


def test_is_excluded_status_live_is_false() -> None:
    """live is not excluded from comparison views."""
    assert is_excluded_status("live") is False


def test_is_excluded_status_completed_is_false() -> None:
    """completed is not excluded from comparison views - only its missing
    metrics are flagged as a legacy gap, the row itself still appears."""
    assert is_excluded_status("completed") is False
