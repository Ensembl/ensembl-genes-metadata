from fastapi import APIRouter, HTTPException

from metadata_app.backend.app.models.gca_lookup_schemas import (
    GcaLookupRequest,
)
from metadata_app.backend.app.services.gca_lookup_service import (
    get_gca_info,
)

gca_lookup = APIRouter()


@gca_lookup.post("/gca_lookup")
def lookup_gca(filters: GcaLookupRequest):
    try:
        return get_gca_info(filters.gca)

    except HTTPException as e:
        raise e

    except Exception:
        import logging

        logging.exception("Unhandled exception in lookup_gca")

        raise HTTPException(
            status_code=500,
            detail="Internal server error",
        )
