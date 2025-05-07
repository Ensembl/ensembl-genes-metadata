from typing import List, Optional, Dict
from pydantic import BaseModel
from datetime import date


class AnnotationFilterRequest(BaseModel):
    bioproject_id: Optional[List[str]] = None
    annotation_date: Optional[date] = None
    taxon_id: Optional[int] = None
    release_type: Optional[List[str]] = None
