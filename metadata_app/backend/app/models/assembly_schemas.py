from typing import List, Optional, Dict
from pydantic import BaseModel
from datetime import date


class AssemblyFilterRequest(BaseModel):
    bioproject_id: Optional[List[str]] = None
    metric_thresholds: Optional[Dict[str, float]] = None
    asm_level: Optional[List[str]] = None
    asm_type: Optional[List[str]] = None
    release_date: Optional[date] = None
    taxon_id: Optional[int] = None
    current: Optional[bool] = False
    pipeline: Optional[List[str]] = None
    transc: Optional[bool] = False