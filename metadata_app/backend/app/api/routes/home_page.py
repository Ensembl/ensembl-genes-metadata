# app/api/routes/home_page.py
import logging
from datetime import date
from fastapi import APIRouter, HTTPException, Depends
from pydantic import BaseModel
from typing import List, Optional

from metadata_app.backend.app.services.home_page_service import get_annotation_counts_by_bioproject
from metadata_app.backend.app.services.home_page_service import get_assemblies_per_year
from metadata_app.backend.app.services.home_page_service import get_annotations_per_year
from metadata_app.backend.app.services.home_page_service import get_metadata_registry_update_dates
from metadata_app.backend.app.services.home_page_service import get_transcriptomic_registry_update_dates
from metadata_app.backend.app.services.home_page_service import get_annotations_per_genebuilder
from metadata_app.backend.app.services.home_page_service import bin_by_genebuild_method




home_page = APIRouter()


class FilterInput(BaseModel):
    bioproject_id: Optional[List[str]] = None
    release_date: Optional[date] = None
    taxon_id: Optional[int] = None
    release_type: Optional[List[str]] = None


@home_page.get("/home/bioproject")
def get_annotation_counts():
    """
    Get counts of annotations per BioProject, with friendly project names.
    """
    return get_annotation_counts_by_bioproject()

@home_page.get("/home/assemblies")
def get_assemblies_year():
    """
    Get counts of assemblies per year.
    """
    return get_assemblies_per_year()

@home_page.get("/home/annotations")
def get_annotations_year():
    """
    Get counts of assemblies per year.
    """
    return get_annotations_per_year()

@home_page.get("/home/meta_update")
def metadata_registry_update_dates():
    """
    Get counts of assemblies per year.
    """
    return get_metadata_registry_update_dates()

@home_page.get("/home/transc_update")
def transcriptomic_registry_update_dates():
    """
    Get counts of assemblies per year.
    """
    return get_transcriptomic_registry_update_dates()

@home_page.get("/home/genebuilder")
def annotations_per_genebuilder():
    """
    Get counts of assemblies per year.
    """
    return get_annotations_per_genebuilder()

@home_page.post("/home/method_summary")
def filter_genebuild_summary(filters: FilterInput):
    try:
        summary_df = bin_by_genebuild_method(
            bioproject_id=filters.bioproject_id,
            release_date=filters.release_date,
            taxon_id=filters.taxon_id,
            release_type=filters.release_type
        )
        return summary_df

    except HTTPException as e:
        raise e  # Preserve FastAPI-specific HTTP exceptions

    except Exception as e:
        logging.exception("Unhandled exception in filter_genebuild_summary")
        raise HTTPException(status_code=500, detail="Internal server error")