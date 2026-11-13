from typing import List, Optional, Dict
from pydantic import BaseModel
from datetime import date
from metadata_app.backend.app.models.filter_schemas import ScopeFilterRequest


class ReportFilterRequest(ScopeFilterRequest):
    start_date: Optional[date] = None
    end_date: Optional[date] = None
