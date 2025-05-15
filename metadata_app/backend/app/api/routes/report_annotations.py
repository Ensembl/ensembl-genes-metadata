from fastapi import APIRouter, HTTPException
from metadata_app.backend.app.models.report_schemas import ReportFilterRequest
from metadata_app.backend.app.services.report_annotation_service import generate_report

report = APIRouter()

@report.post("/report/filter")
def filter_annotations(filters: ReportFilterRequest):
    try:
        anno_wide, number_of_annotations, method_report, num_unique_taxa, top_3_taxa, project_report, average_busco, main_report = generate_report(
            bioproject_id=filters.bioproject_id,
            end_date=filters.end_date,
            start_date=filters.start_date,
            group_name=filters.group_name,
            taxon_id=filters.taxon_id,
            release_type=filters.release_type
        )

        return {
            "main_report": main_report.to_dict(orient="records"),
            "number_of_annotations": number_of_annotations.to_dict(orient="records"),
            "method_report": method_report.to_dict(orient="records"),
            "num_unique_taxa": {"value": num_unique_taxa},
            "top_3_taxa": top_3_taxa.to_dict(orient="records"),
            "project_report": project_report.to_dict(orient="records"),
            "average_busco": {"value": average_busco},

        "downloadables_anno": {
                "anno_main": main_report.to_csv(index=False),
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