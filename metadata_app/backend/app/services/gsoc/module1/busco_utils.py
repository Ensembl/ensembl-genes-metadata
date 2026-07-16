"""
busco_utils.py
Shared utility functions for parsing and interpreting BUSCO score strings.
BUSCO scores in the database are stored as strings like:
"C:94.3%[S:91.2%,D:3.1%],F:2.1%,M:3.6%,n:255"
Newer BUSCO runs may include additional trailing fields after ``n:``, e.g.
"C:97.8%[S:93.6%,D:4.2%],F:0.3%,M:1.9%,n:7207,E:4.5%". These are captured
generically in the ``extra`` dict of the parse result rather than assumed.
This module provides a single canonical parsing location used by both
Module 1 (per-genome reports) and Module 2 (comparative analysis).
"""

import re
from typing import Dict, Optional, Union

from metadata_app.backend.app.services.gsoc.module1.logging_utils import (  # pylint: disable=import-error
    get_logger,
)

logger = get_logger(__name__)

# Canonical quality label vocabulary. This is the single source of truth
# for BUSCO quality bands - all other modules (report_renderer, html_renderer)
# must source their labels/colours from this list rather than redefining
# their own threshold checks or label strings.
QUALITY_THRESHOLDS = (
    ("Excellent", 95.0),
    ("Good", 85.0),
    ("Moderate", 70.0),
)
QUALITY_UNKNOWN = "Unknown"
QUALITY_POOR = "Poor"

# Fields always present in a standard BUSCO short summary string.
_CORE_FIELD_PATTERNS = {
    "complete": r"C:(\d+\.?\d*)%",
    "single_copy": r"S:(\d+\.?\d*)%",
    "duplicated": r"D:(\d+\.?\d*)%",
    "fragmented": r"F:(\d+\.?\d*)%",
    "missing": r"M:(\d+\.?\d*)%",
}
_N_GENES_PATTERN = r"n:(\d+)"

# Any other "Letter:value%" token is captured generically into `extra`
# rather than silently dropped. This covers fields such as "E:4.5%" that
# are not always present and whose exact meaning may vary by BUSCO version.
_GENERIC_EXTRA_PATTERN = re.compile(r"(?<![A-Za-z])([A-Za-z]):(\d+\.?\d*)%")
_KNOWN_LETTERS = {"C", "S", "D", "F", "M"}


def parse_busco_string(
    busco_string: str,
) -> Dict[str, Optional[Union[float, int, Dict[str, float]]]]:
    """
    Parse a BUSCO score string into its component values.

    Args:
        busco_string: Raw BUSCO string e.g. "C:94.3%[S:91.2%,D:3.1%],F:2.1%,M:3.6%,n:255"
            Leading/trailing whitespace and letter case are normalised before
            matching, so "  c:94.3%[s:91.2%...] " is also accepted.

    Returns:
        Dictionary with keys: complete, single_copy, duplicated, fragmented,
        missing, n_genes, extra. All values are floats except n_genes (int)
        and extra (a dict of any additional "Letter:value%" fields found,
        e.g. {"E": 4.5} for "...,n:7207,E:4.5%"). Returns None for fields
        that could not be parsed; extra is always a dict (empty if no
        additional fields were present).

    Example:
        >>> parse_busco_string("C:94.3%[S:91.2%,D:3.1%],F:2.1%,M:3.6%,n:255")
        {'complete': 94.3, 'single_copy': 91.2, 'duplicated': 3.1,
         'fragmented': 2.1, 'missing': 3.6, 'n_genes': 255, 'extra': {}}
    """
    result: Dict[str, Optional[Union[float, int, Dict[str, float]]]] = {
        "complete": None,
        "single_copy": None,
        "duplicated": None,
        "fragmented": None,
        "missing": None,
        "n_genes": None,
        "extra": {},
    }

    if not isinstance(busco_string, str):
        logger.warning("Invalid BUSCO string received: %r", busco_string)
        return result

    normalised = busco_string.strip()
    if not normalised:
        return result

    if not re.search(r"[Cc]:", normalised):
        logger.warning("Invalid BUSCO string received: %r", busco_string)
        return result

    try:
        for key, pattern in _CORE_FIELD_PATTERNS.items():
            match = re.search(pattern, normalised, re.IGNORECASE)
            result[key] = float(match.group(1)) if match else None

        n_genes_match = re.search(_N_GENES_PATTERN, normalised, re.IGNORECASE)
        result["n_genes"] = int(n_genes_match.group(1)) if n_genes_match else None

        extra: Dict[str, float] = {}
        for letter, value in _GENERIC_EXTRA_PATTERN.findall(normalised):
            if letter.upper() in _KNOWN_LETTERS:
                continue
            extra[letter.upper()] = float(value)
        result["extra"] = extra

    except re.error as exc:
        logger.error("Failed to parse BUSCO string '%s': %s", busco_string, exc)

    return result


def get_busco_complete(busco_string: str) -> Optional[float]:
    """
    Convenience function - extract just the complete BUSCO percentage.
    This is the most commonly used single value from the BUSCO string.

    Args:
        busco_string: Raw BUSCO string

    Returns:
        Complete BUSCO percentage as float, or None if unparseable
    """
    value = parse_busco_string(busco_string)["complete"]
    if isinstance(value, (float, int)):
        return float(value)
    return None


def busco_quality_label(complete_pct: Optional[float]) -> str:
    """
    Convert a BUSCO complete percentage to a human-readable quality label.

    This is the single canonical source of BUSCO quality labels. Any code
    that needs to colour-code or style a quality band (e.g. html_renderer,
    report_renderer) should call this function rather than re-implementing
    its own threshold checks, so the label vocabulary never drifts out of
    sync across the codebase.

    Args:
        complete_pct: Complete BUSCO percentage (0-100)

    Returns:
        Quality label string: 'Excellent', 'Good', 'Moderate', 'Poor', or 'Unknown'
    """
    if complete_pct is None:
        return QUALITY_UNKNOWN
    for label, threshold in QUALITY_THRESHOLDS:
        if complete_pct >= threshold:
            return label
    return QUALITY_POOR


# Thresholds for protein vs assembly BUSCO difference flagging.
# Based on Anna's rules (2026-07-08): if protein BUSCO is much lower
# than assembly BUSCO, it suggests the annotation is underperforming
# relative to what the assembly quality would predict.
BUSCO_DIFF_EXCELLENT = 2.0
BUSCO_DIFF_WARNING = 5.0
BUSCO_DIFF_CRITICAL = 10.0

DIFF_LABEL_EXCELLENT = "Excellent"
DIFF_LABEL_ACCEPTABLE = "Acceptable"
DIFF_LABEL_WARNING = "Investigate"
DIFF_LABEL_CRITICAL = "Problematic"
DIFF_LABEL_UNKNOWN = "Unknown"


def busco_diff_label(
    protein_pct: Optional[float],
    assembly_pct: Optional[float],
) -> str:
    """
    Classify the difference between protein and assembly BUSCO completeness.

    When assembly BUSCO is high but protein BUSCO is significantly lower,
    it suggests the annotation may be missing genes that the assembly
    quality would predict should be annotatable.

    Args:
        protein_pct: Protein BUSCO completeness percentage (0-100).
        assembly_pct: Assembly BUSCO completeness percentage (0-100).

    Returns:
        A label string reflecting how far protein BUSCO falls below
        assembly BUSCO. Returns DIFF_LABEL_UNKNOWN if either value
        is None.

    Examples:
        >>> busco_diff_label(97.0, 98.0)
        'Excellent'
        >>> busco_diff_label(92.0, 98.0)
        'Investigate'
        >>> busco_diff_label(85.0, 98.0)
        'Problematic'
    """
    if protein_pct is None or assembly_pct is None:
        return DIFF_LABEL_UNKNOWN

    diff = assembly_pct - protein_pct

    if diff < BUSCO_DIFF_EXCELLENT:
        return DIFF_LABEL_EXCELLENT
    if diff < BUSCO_DIFF_WARNING:
        return DIFF_LABEL_ACCEPTABLE
    if diff < BUSCO_DIFF_CRITICAL:
        return DIFF_LABEL_WARNING
    return DIFF_LABEL_CRITICAL
