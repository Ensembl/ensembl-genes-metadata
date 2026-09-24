# metadata_app/backend/app/api/routes/transcriptomics.py

import logging
import pandas as pd
from typing import Dict, List
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from metadata_app.backend.app.services.transcriptomics_service import (
    get_transcriptomic_assessment_for_ids,
    add_transc_data_to_df,
    get_metadata_from_registry,
)

# Create router for transcriptomics endpoints
transcriptomics = APIRouter(
    tags=["transcriptomics"], responses={404: {"description": "Not found"}}
)


class TranscriptomicRegistryRequest(BaseModel):
    taxon_ids: List[int]


@transcriptomics.post("/assessment")
async def get_transcriptomic_assessment_route(taxonomy_dict: Dict):
    """
    Get transcriptomic assessment data for the provided taxonomy dictionary.

    Args:
            taxonomy_dict: Dictionary mapping lowest taxon IDs to their taxonomy hierarchies

    Returns:
            Dictionary with transcriptomic assessment data
    """
    try:
        trans_df = get_transcriptomic_assessment_for_ids(taxonomy_dict)

        # Convert DataFrame to dictionary for response
        results = []
        for _, row in trans_df.iterrows():
            results.append(
                {
                    "taxon_id": row["taxon_id"],
                    "assessment_date": (
                        row["transc_assess_date"].isoformat()
                        if pd.notnull(row["transc_assess_date"])
                        else None
                    ),
                    "status": row["transc_status"],
                }
            )

        return {"transcriptomic_assessments": results}

    except Exception as e:
        logging.error(f"Error in transcriptomic assessment endpoint: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Internal server error: {str(e)}")


@transcriptomics.post("/enrich")
async def enrich_with_transcriptomic_data(data: Dict):
    """
    Enrich the provided data with transcriptomic assessment information.

    Args:
            data: Dictionary containing info_df (as records) and taxonomy_dict

    Returns:
            Dictionary with enriched data
    """
    try:
        if "info_records" not in data or "taxonomy_dict" not in data:
            raise HTTPException(
                status_code=400,
                detail="Request must include 'info_records' and 'taxonomy_dict'",
            )

        # Convert records to DataFrame
        info_df = pd.DataFrame(data["info_records"])
        taxonomy_dict = data["taxonomy_dict"]

        # Add transcriptomic data
        enriched_df = add_transc_data_to_df(info_df, taxonomy_dict)

        # Convert back to records for response
        return {"enriched_data": enriched_df.to_dict(orient="records")}

    except Exception as e:
        logging.error(f"Error in transcriptomic enrichment endpoint: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Internal server error: {str(e)}")


@transcriptomics.post("/registry")
async def get_transcriptomic_registry_metadata(req: TranscriptomicRegistryRequest):
    try:
        df = get_metadata_from_registry(req.taxon_ids)

        if df.empty:
            return {"records": []}

        for column in ["transc_assess_date"]:
            if column in df.columns:
                df[column] = pd.to_datetime(df[column], errors="coerce").dt.strftime(
                    "%Y-%m-%d"
                )
                df[column] = df[column].where(df[column].notna(), None)

        return {"records": df.where(pd.notna(df), None).to_dict(orient="records")}

    except Exception as e:
        logging.error(f"Error in transcriptomic registry endpoint: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Internal server error: {str(e)}")
