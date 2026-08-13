import logging

import pandas as pd

from metadata_app.backend.app.core.database import get_db_connection


PROJECT_COLUMN_NAMES = [
    "gca",
    "lowest_taxon_id",
    "scientific_name",
    "asm_name",
    "asm_level",
    "gb_status",
    "genebuilder",
    "infra_name",
]


PROJECT_FILTERS = {
    "dtol": {"field": "mb.bioproject_name", "value": "DToL"},
    "tol": {"field": "mb.bioproject_name", "value": "ToL"},
    "erga": {"field": "mb.bioproject_name", "value": "ERGA"},
    "asg": {"field": "mb.bioproject_name", "value": "ASG"},
    "erga-pilot": {"field": "mb.bioproject_name", "value": "ERGA_pilot"},
    "erga-bge": {"field": "mb.bioproject_name", "value": "ERGA/BGE"},
    "vgp": {"field": "mb.bioproject_name", "value": "VGP"},
    "ebp": {"field": "mb.bioproject_name", "value": "EBP"},
    "hprc": {"field": "mb.bioproject_name", "value": "HPRC"},
    "cbp": {"field": "mb.bioproject_name", "value": "CBP"},
    "aegis": {"field": "mb.bioproject_name", "value": "AEGIS"},
    "rodent2k": {"field": "mb.bioproject_name", "value": "Rodent2K"},
    "atlasea": {"field": "mb.bioproject_name", "value": "ATLASea"},
    "laca": {"field": "g.group_name", "value": "LACA"},
}


def _get_project_records(filter_field: str, filter_value: str):
    query = f"""
        SELECT
            CONCAT(a.gca_chain, '.', a.gca_version) AS gca,
            a.lowest_taxon_id,
            s.scientific_name,
            a.asm_name,
            a.asm_level,
            gb.gb_status,
            gb.genebuilder,
            o.infra_name
        FROM assembly a
        LEFT JOIN genebuild_status gb ON gb.assembly_id = a.assembly_id
        LEFT JOIN bioproject b ON a.assembly_id = b.assembly_id
        LEFT JOIN main_bioproject mb ON b.bioproject_id = mb.bioproject_id
        LEFT JOIN species s ON a.lowest_taxon_id = s.lowest_taxon_id
        LEFT JOIN organism o ON a.assembly_id = o.assembly_id
        LEFT JOIN custom_group g
            ON (
                (g.group_type = 'taxon' AND a.lowest_taxon_id = g.item)
                OR
                (g.group_type = 'assembly' AND a.gca_chain = g.item)
            )
        WHERE {filter_field} = %s
    """

    with get_db_connection("meta") as conn:
        cursor = conn.cursor()
        cursor.execute(query, (filter_value,))
        return cursor.fetchall()


def _format_project_result(result):
    df = pd.DataFrame(result, columns=PROJECT_COLUMN_NAMES)
    if df.empty:
        return []

    df["gb_status"] = df["gb_status"].fillna("not_started")
    df = df.drop(columns=["asm_name"])
    df = df.drop_duplicates(subset=["gca", "gb_status"], keep="first")
    return df.to_dict(orient="records")


def get_project(project_key: str):
    config = PROJECT_FILTERS.get(project_key)
    if not config:
        logging.error("Unknown project key: %s", project_key)
        return []

    try:
        result = _get_project_records(
            filter_field=config["field"],
            filter_value=config["value"],
        )
        return _format_project_result(result)
    except Exception as exc:
        logging.error("Error fetching project records for %s: %s", project_key, exc)
        return []


def get_dtol():
    return get_project("dtol")


def get_tol():
    return get_project("tol")


def get_erga():
    return get_project("erga")


def get_asg():
    return get_project("asg")


def get_erga_pilot():
    return get_project("erga-pilot")


def get_erga_bge():
    return get_project("erga-bge")


def get_vgp():
    return get_project("vgp")


def get_ebp():
    return get_project("ebp")


def get_hprc():
    return get_project("hprc")


def get_cbp():
    return get_project("cbp")


def get_laca():
    return get_project("laca")


def get_aegis():
    return get_project("aegis")


def get_rodent2k():
    return get_project("rodent2k")


def get_atlasea():
    return get_project("atlasea")
