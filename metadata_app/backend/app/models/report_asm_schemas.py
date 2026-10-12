from typing import List, Optional, Dict
from pydantic import BaseModel
from datetime import date


class ReportFilterRequest(BaseModel):
    bioproject_id: Optional[List[str]] = None
    group_name: Optional[List[str]] = None
    start_date: Optional[date] = None
    end_date: Optional[date] = None
    taxon_id: Optional[List[int]] = None
    candidate: Optional[bool] = False
    transc: Optional[bool] = False
    transc_ena: Optional[bool] = False
    non_annotated: Optional[bool] = True
