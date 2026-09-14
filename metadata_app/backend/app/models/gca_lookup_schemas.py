from pydantic import BaseModel
from typing import Optional


class GcaLookupRequest(BaseModel):
    gca: list[str]
