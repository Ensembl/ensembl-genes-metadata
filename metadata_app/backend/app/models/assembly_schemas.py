from typing import Optional
from pydantic import BaseModel
from datetime import date
from metadata_app.backend.app.models.filter_schemas import (
    AssemblyDimensionFilters,
    ScopeFilterRequest,
)


class AssemblyFilterRequest(ScopeFilterRequest, AssemblyDimensionFilters):
    release_date: Optional[date] = None
    current: Optional[bool] = False
    transc: Optional[bool] = False
    transc_ena: Optional[bool] = False
    non_annotated: Optional[bool] = False
    candidate: Optional[bool] = False
    include_reassessment: Optional[bool] = False
