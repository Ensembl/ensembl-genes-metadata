from typing import List, Optional, Dict
from pydantic import BaseModel
from datetime import date


class ReportFilterRequest(BaseModel):
    bioproject_id: Optional[List[str]] = None
    group_name: Optional[str] = None
    end_date: Optional[date] = None
    start_date: Optional[date] = None
    taxon_id: Optional[List[int]] = None
    release_type: Optional[str] = None
