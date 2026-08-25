# app/api/routes/db_clean.py
import io
import logging
from fastapi import APIRouter, HTTPException, Form
from fastapi.responses import StreamingResponse

from metadata_app.backend.app.services.db_clean_service import server_clean_main

db_clean_router = APIRouter()

@db_clean_router.post("/db_clean/genebuilder")
def get_server_clean_main(genebuilder: str = Form(...)):
    try:
        _, commands = server_clean_main(genebuilder)
        if not commands:
            raise HTTPException(status_code=404, detail="No cleanup commands returned")

        file_content = io.BytesIO(commands.encode("utf-8"))
        headers = {
            "Content-Disposition": f'attachment; filename="cleanup_{genebuilder}.sql"'
        }

        return StreamingResponse(
            file_content,
            media_type="application/sql",
            headers=headers
        )

    except HTTPException:
        raise
    except Exception:
        logging.exception("Unhandled exception in server_clean_main")
        raise HTTPException(status_code=500, detail="Internal server error")