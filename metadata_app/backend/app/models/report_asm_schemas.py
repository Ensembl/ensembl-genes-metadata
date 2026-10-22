from typing import Optional
from datetime import date
from metadata_app.backend.app.models.filter_schemas import (
    AssemblyDimensionFilters,
    ScopeFilterRequest,
)


class ReportFilterRequest(ScopeFilterRequest, AssemblyDimensionFilters):
    start_date: Optional[date] = None
    end_date: Optional[date] = None
    candidate: Optional[bool] = False
    include_reassessment: Optional[bool] = False
    current: Optional[bool] = False
    transc: Optional[bool] = False
    transc_ena: Optional[bool] = False
    non_annotated: Optional[bool] = False
