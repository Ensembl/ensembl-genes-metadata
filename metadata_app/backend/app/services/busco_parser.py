from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Optional


_FULL_BUSCO_RE = re.compile(
    r"""
    C:(?P<complete>\d+(?:\.\d+)?)%
    (?:\[
        S:(?P<single>\d+(?:\.\d+)?)%,
        D:(?P<duplicated>\d+(?:\.\d+)?)%
    \])?
    (?:,F:(?P<fragmented>\d+(?:\.\d+)?)%)?
    (?:,M:(?P<missing>\d+(?:\.\d+)?)%)?
    (?:,n:(?P<searched>\d+))?
    """,
    re.VERBOSE,
)


@dataclass(frozen=True)
class BuscoScores:
    complete: Optional[float] = None
    single: Optional[float] = None
    duplicated: Optional[float] = None
    fragmented: Optional[float] = None
    missing: Optional[float] = None
    searched: Optional[int] = None


def parse_busco_string(value: object) -> BuscoScores:
    """Parse a BUSCO summary string into typed components.

    Supports both the full string format, e.g.
    ``C:98.2%[S:97.5%,D:0.7%],F:0.8%,M:1.0%,n:255``, and older
    completeness-only values such as ``C:98.2%``.
    """

    if value is None:
        return BuscoScores()

    if not isinstance(value, str):
        value = str(value)

    value = value.strip()
    if not value:
        return BuscoScores()

    match = _FULL_BUSCO_RE.search(value)
    if not match:
        return BuscoScores()

    groups = match.groupdict()
    return BuscoScores(
        complete=_to_float(groups["complete"]),
        single=_to_float(groups["single"]),
        duplicated=_to_float(groups["duplicated"]),
        fragmented=_to_float(groups["fragmented"]),
        missing=_to_float(groups["missing"]),
        searched=_to_int(groups["searched"]),
    )


def _to_float(value: Optional[str]) -> Optional[float]:
    if value is None:
        return None
    return float(value)


def _to_int(value: Optional[str]) -> Optional[int]:
    if value is None:
        return None
    return int(value)
