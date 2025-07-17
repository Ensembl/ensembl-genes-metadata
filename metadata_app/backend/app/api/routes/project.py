# app/api/routes/home_page.py

from fastapi import APIRouter


from metadata_app.backend.app.services.project_service import get_erga
from metadata_app.backend.app.services.project_service import get_vgp
from metadata_app.backend.app.services.project_service import get_ebp
from metadata_app.backend.app.services.project_service import get_erga_bge
from metadata_app.backend.app.services.project_service import get_asg
from metadata_app.backend.app.services.project_service import get_cbp
from metadata_app.backend.app.services.project_service import get_dtol
from metadata_app.backend.app.services.project_service import get_hprc
from metadata_app.backend.app.services.project_service import get_erga_pilot


project_router = APIRouter()


@project_router.get("/project/erga")
def get_erga_project():
    """
    """
    return get_erga()

@project_router.get("/project/vgp")
def get_vgp_project():
    return get_vgp()

@project_router.get("/project/ebp")
def get_ebp_project():
    return get_ebp()

@project_router.get("/project/erga-bge")
def get_erga_bge_project():
    return get_erga_bge()

@project_router.get("/project/asg")
def get_asg_project():
    return get_asg()

@project_router.get("/project/cbp")
def get_cbp_project():
    return get_cbp()

@project_router.get("/project/dtol")
def get_dtol_project():
    return get_dtol()

@project_router.get("/project/hprc")
def get_hprc_project():
    return get_hprc()

@project_router.get("/project/erga-pilot")
def get_erga_pilot_project():
    return get_erga_pilot()



