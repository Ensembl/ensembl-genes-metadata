# app/api/routes/handover.py
import logging
from datetime import date
from fastapi import APIRouter, HTTPException, Depends
from pydantic import BaseModel
from typing import List, Optional

from metadata_app.backend.app.services.handover_service import get_ready_to_ho

handover_router = APIRouter()

@handover_router.get("/handover/genebuilder")
def get_ready_to_handover():
    """
    Get HO ready GCAs per genebuilder.
    """
    return get_ready_to_ho()
