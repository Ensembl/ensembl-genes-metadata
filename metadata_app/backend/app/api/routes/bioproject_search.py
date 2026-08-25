from fastapi import APIRouter, HTTPException, Query
from metadata_app.backend.app.services.bioproject_search_service import search_bioproject

router = APIRouter()

@router.get("/bioprojects/search")
def bioproject_search(q: str = Query(..., description="Text to search by name or ID")):
    return search_bioproject(q)


