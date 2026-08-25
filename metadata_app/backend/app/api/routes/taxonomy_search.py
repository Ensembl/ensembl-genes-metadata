from fastapi import APIRouter, HTTPException, Query
from metadata_app.backend.app.services.taxonomy_search_service import search_taxonomy

router = APIRouter()

@router.get("/taxonomy/search")
def taxonomy_search(q: str = Query(..., description="Text to search by taxonomy name or id")):
    return search_taxonomy(q)


