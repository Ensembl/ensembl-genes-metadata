# app/api/routes/handover.py
import logging
from datetime import date
from fastapi import APIRouter, HTTPException, Depends
from pydantic import BaseModel
from typing import List, Optional

from metadata_app.backend.app.services.handover_service import (
    get_ready_to_ho,
    update_gca,
)

handover_router = APIRouter()


class GenebuilderRequest(BaseModel):
    genebuilder: str


class HandoverItem(BaseModel):
    gca: str
    annotation_method: str


class UpdateGCARequest(BaseModel):
    genebuilder: str
    items: List[HandoverItem]
    new_status: str


@handover_router.post("/handover/genebuilder")
async def get_ready_to_handover(req: GenebuilderRequest):
    genebuilder = req.genebuilder
    result = get_ready_to_ho(genebuilder)

    if not result:
        raise HTTPException(status_code=404, detail="No results found")

    if isinstance(result, tuple) and isinstance(result[0], str):
        raise HTTPException(status_code=400, detail=result[0])

    (
        df_ready,
        count_ho_ready,
        count_data,
        count_pending,
        list_data,
        list_pending,
        annotation_overview,
        status_summary,
    ) = result

    return {
        "df_ready": df_ready,
        "count_ho_ready": count_ho_ready,
        "count_data": count_data,
        "count_pending": count_pending,
        "list_data": list_data,
        "list_pending": list_pending,
        "annotation_overview": annotation_overview,
        "status_summary": status_summary,
    }


@handover_router.post("/handover/change_status")
async def abandon_gcas(req: UpdateGCARequest):
    if not req.items:
        raise HTTPException(status_code=400, detail="No GCAs provided")

    updated_rows = update_gca(
        genebuilder=req.genebuilder,
        items=[item.dict() for item in req.items],
        new_status=req.new_status,
    )

    if updated_rows == 0:
        raise HTTPException(
            status_code=404,
            detail="No records were updated (check GCA, annotation_method, or status)",
        )

    return {
        "status": "ok",
        "updated_rows": updated_rows,
    }
