from typing import List, Optional, Dict
from pydantic import BaseModel
from datetime import date
from metadata_app.backend.app.models.filter_schemas import ScopeFilterRequest


class AnnotationFilterRequest(ScopeFilterRequest):
    annotation_date: Optional[date] = None
    bioproject_name: Optional[List[str]] = None
    gca: Optional[List[str]] = None
