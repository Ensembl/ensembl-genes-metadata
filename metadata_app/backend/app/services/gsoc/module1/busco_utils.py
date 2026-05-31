"""
busco_utils.py
Shared utility functions for parsing and interpreting BUSCO score strings.
BUSCO scores in the database are stored as strings like:
"C:94.3%[S:91.2%,D:3.1%],F:2.1%,M:3.6%,n:255"
This module provides a single canonical parsing location used by both
Module 1 (per-genome reports) and Module 2 (comparative analysis).
"""

import re
import logging
from typing import Optional


def parse_busco_string(busco_string: str) -> dict[str, Optional[float | int]]:
    """
    Parse a BUSCO score string into its component values.

    Args:
        busco_string: Raw BUSCO string e.g. "C:94.3%[S:91.2%,D:3.1%],F:2.1%,M:3.6%,n:255"

    Returns:
        Dictionary with keys: complete, single_copy, duplicated, fragmented, missing, n_genes
        All float values except n_genes (int). Returns None for unparseable fields.

    Example:
        >>> parse_busco_string("C:94.3%[S:91.2%,D:3.1%],F:2.1%,M:3.6%,n:255")
        {'complete': 94.3, 'single_copy': 91.2, 'duplicated': 3.1,
         'fragmented': 2.1, 'missing': 3.6, 'n_genes': 255}
    """
    result: dict[str, Optional[float | int]] = {
        "complete": None,
        "single_copy": None,
        "duplicated": None,
        "fragmented": None,
        "missing": None,
        "n_genes": None,
    }

    if not busco_string or not isinstance(busco_string, str):
        logging.warning("Invalid BUSCO string received: %s", busco_string)
        return result

    try:
        complete_match = re.search(r"C:(\d+\.?\d*)%", busco_string)
        single_match = re.search(r"S:(\d+\.?\d*)%", busco_string)
        duplicated_match = re.search(r"D:(\d+\.?\d*)%", busco_string)
        fragmented_match = re.search(r"F:(\d+\.?\d*)%", busco_string)
        missing_match = re.search(r"M:(\d+\.?\d*)%", busco_string)
        n_genes_match = re.search(r"n:(\d+)", busco_string)

        result["complete"] = float(complete_match.group(1)) if complete_match else None
        result["single_copy"] = float(single_match.group(1)) if single_match else None
        result["duplicated"] = (
            float(duplicated_match.group(1)) if duplicated_match else None
        )
        result["fragmented"] = (
            float(fragmented_match.group(1)) if fragmented_match else None
        )
        result["missing"] = float(missing_match.group(1)) if missing_match else None
        result["n_genes"] = int(n_genes_match.group(1)) if n_genes_match else None

    except re.error as e:
        logging.error("Failed to parse BUSCO string '%s': %s", busco_string, e)

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
    return float(value) if value is not None else None


def busco_quality_label(complete_pct: Optional[float]) -> str:
    """
    Convert a BUSCO complete percentage to a human-readable quality label.

    Args:
        complete_pct: Complete BUSCO percentage (0-100)

    Returns:
        Quality label string: 'Excellent', 'Good', 'Moderate', 'Poor', or 'Unknown'
    """
    if complete_pct is None:
        return "Unknown"
    if complete_pct >= 95:
        return "Excellent"
    if complete_pct >= 85:
        return "Good"
    if complete_pct >= 70:
        return "Moderate"
    return "Poor"
