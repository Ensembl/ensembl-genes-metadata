from typing import Dict, List, Optional

from pydantic import BaseModel


class ScopeFilterRequest(BaseModel):
    """Filters shared by assembly and annotation endpoints."""

    bioproject_id: Optional[List[str]] = None
    group_name: Optional[List[str]] = None
    taxon_id: Optional[List[int]] = None


class AssemblyDimensionFilters(BaseModel):
    """Optional assembly filters shared by assembly and report APIs."""

    metric_thresholds: Optional[Dict[str, float]] = None
    asm_level: Optional[List[str]] = None
    asm_type: Optional[List[str]] = None
    pipeline: Optional[List[str]] = None
    gca: Optional[List[str]] = None
