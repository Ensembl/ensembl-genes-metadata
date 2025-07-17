from fastapi import APIRouter, HTTPException
from metadata_app.backend.app.models.annotation_schemas import AnnotationFilterRequest
from metadata_app.backend.app.services.annotations_service import generate_tables

annotations = APIRouter()

@annotations.post("/annotations/filter")
def filter_annotations(filters: AnnotationFilterRequest):
    try:
        anno_wide, anno_main = generate_tables(
            bioproject_id=filters.bioproject_id,
            annotation_date=filters.annotation_date,
            taxon_id=filters.taxon_id,
            release_type=filters.release_type,
        )

        return {
            "anno_main": anno_main.to_dict(orient="records"),
            "downloadables_anno": {
                "anno_main": anno_main.to_csv(index=False),
                "anno_wide": anno_wide.to_csv(index=False)
            }
        }

    except HTTPException as e:
        raise e  # re-raise to return proper status like 404
    except Exception as e:
        # Log the actual exception
        import logging
        logging.exception("Unhandled exception in filter_annotations")
        raise HTTPException(status_code=500, detail="Internal server error")