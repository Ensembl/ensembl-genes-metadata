import logging
import json
from pathlib import Path

from fastapi import APIRouter, HTTPException

from metadata_app.backend.app.models.qc_schemas import QcFilterRequest
from metadata_app.backend.app.services.qc_service import cohort_feature_means, prepare_qc_records, query_qc_data

qc = APIRouter()
CLADE_SETTINGS_PATH = Path(__file__).resolve().parents[3] / "data" / "clade_settings.json"


@qc.get("/clades")
def list_qc_clades():
    """Return the clade names configured for annotation QC."""
    try:
        with CLADE_SETTINGS_PATH.open("r", encoding="utf-8") as handle:
            clades = sorted(json.load(handle).keys(), key=str.casefold)
        return {"clades": clades}
    except Exception as exc:
        logging.exception("Unable to load QC clade settings")
        raise HTTPException(status_code=500, detail="Unable to load clade settings") from exc


@qc.post("/annotation/filter")
def filter_annotation_qc(filters: QcFilterRequest):
    """Return live annotation-QC records for the requested cohort."""
    try:
        qc_df = query_qc_data(
            annotation_date=filters.annotation_date,
            taxon_id=filters.taxon_id,
            bioproject_id=filters.bioproject_id,
            group_name=filters.group_name,
            gca=filters.gca,
            clade=filters.clade,
            status=filters.status,
            annotation_method=filters.annotation_method,
        )
        prepared_records = prepare_qc_records(qc_df)
        cohort_means = cohort_feature_means(prepared_records, filters.gca)
        records = prepared_records
        if filters.gca:
            requested_gcas = {str(value).strip().casefold() for value in filters.gca}
            records = [
                record
                for record in records
                if str(record.get("gca", "")).casefold() in requested_gcas
            ]
        if filters.status:
            records = [record for record in records if record.get("status") == filters.status]
        if filters.clade:
            requested_clade = filters.clade.casefold()
            records = [
                record
                for record in records
                if str(record.get("taxonomy", "")).casefold() == requested_clade
            ]
        return {
            "records": records,
            "cohort_records": prepared_records,
            "count": len(records),
            "statuses": sorted({record["status"] for record in records if record.get("status")}),
            "taxonomies": sorted({record["taxonomy"] for record in records if record.get("taxonomy")}),
            "cohort_means": cohort_means,
        }
    except HTTPException:
        raise
    except Exception as exc:
        logging.exception("Unhandled exception in annotation QC filter")
        raise HTTPException(status_code=500, detail="Internal server error") from exc
