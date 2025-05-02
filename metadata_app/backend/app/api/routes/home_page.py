# app/api/routes/home_page.py

from fastapi import APIRouter
from metadata_app.backend.app.services.home_page_service import get_annotation_counts_by_bioproject
from metadata_app.backend.app.services.home_page_service import get_assemblies_per_year
from metadata_app.backend.app.services.home_page_service import get_annotations_per_year
from metadata_app.backend.app.services.home_page_service import get_metadata_registry_update_dates
from metadata_app.backend.app.services.home_page_service import get_transcriptomic_registry_update_dates
from metadata_app.backend.app.services.home_page_service import get_annotations_per_genebuilder



home_page = APIRouter()


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