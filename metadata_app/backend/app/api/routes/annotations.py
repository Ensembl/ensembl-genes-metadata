from fastapi import APIRouter, HTTPException
from metadata_app.backend.app.models.annotation_schemas import AnnotationFilterRequest
from metadata_app.backend.app.services.annotations_service import generate_tables, compute_gene_stats

annotations = APIRouter()

@annotations.post("/annotations/filter")
def filter_annotations(filters: AnnotationFilterRequest):
    try:
        anno_wide, anno_main = generate_tables(
            bioproject_id=filters.bioproject_id,
            annotation_date=filters.annotation_date,
            taxon_id=filters.taxon_id,
            group_name= filters.group_name,
        )
        
        # Calculate summary statistics
        gene_stats = compute_gene_stats(anno_wide)

        # Build an informative suffix
        suffix = "_".join(
            filter(None, [
                f"bioproject_{'_'.join(filters.bioproject_id)}" if filters.bioproject_id else None,
                f"taxon_{'_'.join(map(str, filters.taxon_id))}" if filters.taxon_id else None,
                f"date_{filters.annotation_date}" if filters.annotation_date else None,
                f"group_{'_'.join(filters.group_name)}" if filters.group_name else None,
            ])
        )

        return {
            "anno_main": anno_main.to_dict(orient="records"),
            "gene_stats": gene_stats,   # 👈 NEW FEATURE
            "downloadables_anno": {
                "anno_main": {
                    "filename": f"anno_main_{suffix}.csv",
                    "csv": anno_main.to_csv(index=False),
                },
                "anno_wide": {
                    "filename": f"anno_wide_{suffix}.csv",
                    "csv": anno_wide.to_csv(index=False),
                },
            },
        }

    except HTTPException as e:
        raise e  # re-raise to return proper status like 404
    except Exception as e:
        # Log the actual exception
        import logging
        logging.exception("Unhandled exception in filter_annotations")
        raise HTTPException(status_code=500, detail="Internal server error")