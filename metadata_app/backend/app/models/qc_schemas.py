from datetime import date
from typing import List, Optional

from pydantic import BaseModel


class QcFilterRequest(BaseModel):
    bioproject_id: Optional[List[str]] = None
    annotation_date: Optional[date] = None
    taxon_id: Optional[List[int]] = None
    group_name: Optional[List[str]] = None
    gca: Optional[List[str]] = None
    status: Optional[str] = None
    clade: Optional[str] = None
    annotation_method: Optional[List[str]] = None
