# app/services/project_service.py
import logging
from metadata_app.backend.app.core.database import get_db_connection
import pandas as pd

def get_dtol():
    """Returns GCA info"""
    try:
        with get_db_connection("meta") as conn:
            cursor = conn.cursor()
            query = """
                SELECT CONCAT(a.gca_chain, '.', a.gca_version) AS gca, a.lowest_taxon_id, s.scientific_name, a.asm_name,
                        a.asm_level, gb.gb_status, gb.genebuilder
                FROM assembly a
                LEFT JOIN genebuild gb on gb.assembly_id = a.assembly_id
                LEFT JOIN bioproject b on a.assembly_id = b.assembly_id
                LEFT JOIN species s ON a.lowest_taxon_id = s.lowest_taxon_id
                LEFT JOIN group_assembly g ON a.assembly_id = g.assembly_id
                JOIN main_bioproject mb ON b.bioproject_id = mb.bioproject_id
                WHERE mb.bioproject_name = 'DToL'
            """
            cursor.execute(query)
            result = cursor.fetchall()

        df = pd.DataFrame(result, columns=[
            "gca",
            "lowest_taxon_id",
            "scientific_name",
            "asm_name",
            "asm_level",
            "gb_status",
            "genebuilder"
        ])
        df["gb_status"] = df["gb_status"].fillna("not_started")
        df = df[~df["asm_name"].str.contains("alternate", case=False, na=False)]
        df = df.drop(columns=["asm_name"])
        df = df[~df["asm_level"].str.lower().isin(["contig", "scaffold"])]

        return df.to_dict(orient="records")

    except Exception as e:
        logging.error(f"Error fetching annotation counts: {e}")
        return []

def get_erga():
    """Returns GCA info"""
    try:
        with get_db_connection("meta") as conn:
            cursor = conn.cursor()
            query = """
                SELECT CONCAT(a.gca_chain, '.', a.gca_version) AS gca, a.lowest_taxon_id, s.scientific_name, a.asm_name,
                        a.asm_level, gb.gb_status, gb.genebuilder
                FROM genebuild gb
                LEFT JOIN assembly a on gb.assembly_id = a.assembly_id
                LEFT JOIN bioproject b on a.assembly_id = b.assembly_id
                LEFT JOIN species s ON a.lowest_taxon_id = s.lowest_taxon_id
                LEFT JOIN group_assembly g ON a.assembly_id = g.assembly_id
                JOIN main_bioproject mb ON b.bioproject_id = mb.bioproject_id
                WHERE mb.bioproject_name = 'ERGA'
            """
            cursor.execute(query)
            result = cursor.fetchall()

        df = pd.DataFrame(result, columns=[
            "gca",
            "lowest_taxon_id",
            "scientific_name",
            "asm_name",
            "asm_level",
            "gb_status",
            "genebuilder"
        ])
        df["gb_status"] = df["gb_status"].fillna("not_started")
        df = df[~df["asm_name"].str.contains("alternate", case=False, na=False)]
        df = df.drop(columns=["asm_name"])
        df = df[~df["asm_level"].str.lower().isin(["contig", "scaffold"])]
        return df.to_dict(orient="records")

    except Exception as e:
        logging.error(f"Error fetching annotation counts: {e}")
        return []


def get_asg():
    """Returns GCA info"""
    try:
        with get_db_connection("meta") as conn:
            cursor = conn.cursor()
            query = """
                SELECT CONCAT(a.gca_chain, '.', a.gca_version) AS gca, a.lowest_taxon_id, s.scientific_name, a.asm_name,
                        a.asm_level, gb.gb_status, gb.genebuilder
                FROM genebuild gb
                LEFT JOIN assembly a on gb.assembly_id = a.assembly_id
                LEFT JOIN bioproject b on a.assembly_id = b.assembly_id
                LEFT JOIN species s ON a.lowest_taxon_id = s.lowest_taxon_id
                LEFT JOIN group_assembly g ON a.assembly_id = g.assembly_id
                JOIN main_bioproject mb ON b.bioproject_id = mb.bioproject_id
                WHERE mb.bioproject_name = 'ASG'
            """
            cursor.execute(query)
            result = cursor.fetchall()

        df = pd.DataFrame(result, columns=[
            "gca",
            "lowest_taxon_id",
            "scientific_name",
            "asm_name",
            "asm_level",
            "gb_status",
            "genebuilder"
        ])
        df["gb_status"] = df["gb_status"].fillna("not_started")
        df = df[~df["asm_name"].str.contains("alternate", case=False, na=False)]
        df = df.drop(columns=["asm_name"])
        df = df[~df["asm_level"].str.lower().isin(["contig", "scaffold"])]
        return df.to_dict(orient="records")

    except Exception as e:
        logging.error(f"Error fetching annotation counts: {e}")
        return []


def get_erga_pilot():
    """Returns GCA info"""
    try:
        with get_db_connection("meta") as conn:
            cursor = conn.cursor()
            query = """
                SELECT CONCAT(a.gca_chain, '.', a.gca_version) AS gca, a.lowest_taxon_id, s.scientific_name, a.asm_name,
                        a.asm_level, gb.gb_status, gb.genebuilder
                FROM genebuild gb
                LEFT JOIN assembly a on gb.assembly_id = a.assembly_id
                LEFT JOIN bioproject b on a.assembly_id = b.assembly_id
                LEFT JOIN species s ON a.lowest_taxon_id = s.lowest_taxon_id
                LEFT JOIN group_assembly g ON a.assembly_id = g.assembly_id
                JOIN main_bioproject mb ON b.bioproject_id = mb.bioproject_id
                WHERE mb.bioproject_name = 'ERGA_pilot'
            """
            cursor.execute(query)
            result = cursor.fetchall()

        df = pd.DataFrame(result, columns=[
            "gca",
            "lowest_taxon_id",
            "scientific_name",
            "asm_name",
            "asm_level",
            "gb_status",
            "genebuilder"
        ])
        df["gb_status"] = df["gb_status"].fillna("not_started")
        df = df[~df["asm_name"].str.contains("alternate", case=False, na=False)]
        df = df.drop(columns=["asm_name"])
        df = df[~df["asm_level"].str.lower().isin(["contig", "scaffold"])]
        return df.to_dict(orient="records")

    except Exception as e:
        logging.error(f"Error fetching annotation counts: {e}")
        return []

def get_erga_bge():
    """Returns GCA info"""
    try:
        with get_db_connection("meta") as conn:
            cursor = conn.cursor()
            query = """
                SELECT CONCAT(a.gca_chain, '.', a.gca_version) AS gca, a.lowest_taxon_id, s.scientific_name, a.asm_name,
                        a.asm_level, gb.gb_status, gb.genebuilder
                FROM genebuild gb
                LEFT JOIN assembly a on gb.assembly_id = a.assembly_id
                LEFT JOIN bioproject b on a.assembly_id = b.assembly_id
                LEFT JOIN species s ON a.lowest_taxon_id = s.lowest_taxon_id
                LEFT JOIN group_assembly g ON a.assembly_id = g.assembly_id
                JOIN main_bioproject mb ON b.bioproject_id = mb.bioproject_id
                WHERE mb.bioproject_name = 'ERGA/BGE'
            """
            cursor.execute(query)
            result = cursor.fetchall()

        df = pd.DataFrame(result, columns=[
            "gca",
            "lowest_taxon_id",
            "scientific_name",
            "asm_name",
            "asm_level",
            "gb_status",
            "genebuilder"
        ])
        df["gb_status"] = df["gb_status"].fillna("not_started")
        df = df[~df["asm_name"].str.contains("alternate", case=False, na=False)]
        df = df.drop(columns=["asm_name"])
        df = df[~df["asm_level"].str.lower().isin(["contig", "scaffold"])]
        return df.to_dict(orient="records")

    except Exception as e:
        logging.error(f"Error fetching annotation counts: {e}")
        return []

def get_vgp():
    """Returns GCA info"""
    try:
        with get_db_connection("meta") as conn:
            cursor = conn.cursor()
            query = """
                SELECT CONCAT(a.gca_chain, '.', a.gca_version) AS gca, a.lowest_taxon_id, s.scientific_name, a.asm_name,
                        a.asm_level, gb.gb_status, gb.genebuilder
                FROM genebuild gb
                LEFT JOIN assembly a on gb.assembly_id = a.assembly_id
                LEFT JOIN bioproject b on a.assembly_id = b.assembly_id
                LEFT JOIN species s ON a.lowest_taxon_id = s.lowest_taxon_id
                LEFT JOIN group_assembly g ON a.assembly_id = g.assembly_id
                JOIN main_bioproject mb ON b.bioproject_id = mb.bioproject_id
                WHERE mb.bioproject_name = 'VGP'
            """
            cursor.execute(query)
            result = cursor.fetchall()

        df = pd.DataFrame(result, columns=[
            "gca",
            "lowest_taxon_id",
            "scientific_name",
            "asm_name",
            "asm_level",
            "gb_status",
            "genebuilder"
        ])
        df["gb_status"] = df["gb_status"].fillna("not_started")
        df = df[~df["asm_name"].str.contains("alternate", case=False, na=False)]
        df = df.drop(columns=["asm_name"])
        df = df[~df["asm_level"].str.lower().isin(["contig", "scaffold"])]
        return df.to_dict(orient="records")

    except Exception as e:
        logging.error(f"Error fetching annotation counts: {e}")
        return []

def get_ebp():
    """Returns GCA info"""
    try:
        with get_db_connection("meta") as conn:
            cursor = conn.cursor()
            query = """
                SELECT CONCAT(a.gca_chain, '.', a.gca_version) AS gca, a.lowest_taxon_id, s.scientific_name, a.asm_name,
                        a.asm_level, gb.gb_status, gb.genebuilder
                FROM genebuild gb
                LEFT JOIN assembly a on gb.assembly_id = a.assembly_id
                LEFT JOIN bioproject b on a.assembly_id = b.assembly_id
                LEFT JOIN species s ON a.lowest_taxon_id = s.lowest_taxon_id
                LEFT JOIN group_assembly g ON a.assembly_id = g.assembly_id
                JOIN main_bioproject mb ON b.bioproject_id = mb.bioproject_id
                WHERE mb.bioproject_name = 'EBP'
            """
            cursor.execute(query)
            result = cursor.fetchall()

        df = pd.DataFrame(result, columns=[
            "gca",
            "lowest_taxon_id",
            "scientific_name",
            "asm_name",
            "asm_level",
            "gb_status",
            "genebuilder"
        ])
        df["gb_status"] = df["gb_status"].fillna("not_started")
        df = df[~df["asm_name"].str.contains("alternate", case=False, na=False)]
        df = df.drop(columns=["asm_name"])
        df = df[~df["asm_level"].str.lower().isin(["contig", "scaffold"])]
        return df.to_dict(orient="records")

    except Exception as e:
        logging.error(f"Error fetching annotation counts: {e}")
        return []

def get_hprc():
    """Returns GCA info"""
    try:
        with get_db_connection("meta") as conn:
            cursor = conn.cursor()
            query = """
                SELECT CONCAT(a.gca_chain, '.', a.gca_version) AS gca, a.lowest_taxon_id, s.scientific_name, a.asm_name,
                        a.asm_level, gb.gb_status, gb.genebuilder
                FROM genebuild gb
                LEFT JOIN assembly a on gb.assembly_id = a.assembly_id
                LEFT JOIN bioproject b on a.assembly_id = b.assembly_id
                LEFT JOIN species s ON a.lowest_taxon_id = s.lowest_taxon_id
                LEFT JOIN group_assembly g ON a.assembly_id = g.assembly_id
                JOIN main_bioproject mb ON b.bioproject_id = mb.bioproject_id
                WHERE mb.bioproject_name = 'HPRC'
            """
            cursor.execute(query)
            result = cursor.fetchall()

        df = pd.DataFrame(result, columns=[
            "gca",
            "lowest_taxon_id",
            "scientific_name",
            "asm_name",
            "asm_level",
            "gb_status",
            "genebuilder"
        ])
        df["gb_status"] = df["gb_status"].fillna("not_started")
        df = df[~df["asm_name"].str.contains("alternate", case=False, na=False)]
        df = df.drop(columns=["asm_name"])
        df = df[~df["asm_level"].str.lower().isin(["contig", "scaffold"])]
        return df.to_dict(orient="records")

    except Exception as e:
        logging.error(f"Error fetching annotation counts: {e}")
        return []

def get_cbp():
    """Returns GCA info"""
    try:
        with get_db_connection("meta") as conn:
            cursor = conn.cursor()
            query = """
                SELECT CONCAT(a.gca_chain, '.', a.gca_version) AS gca, a.lowest_taxon_id, s.scientific_name, a.asm_name,
                        a.asm_level, gb.gb_status, gb.genebuilder
                FROM genebuild gb
                LEFT JOIN assembly a on gb.assembly_id = a.assembly_id
                LEFT JOIN bioproject b on a.assembly_id = b.assembly_id
                LEFT JOIN species s ON a.lowest_taxon_id = s.lowest_taxon_id
                LEFT JOIN group_assembly g ON a.assembly_id = g.assembly_id
                JOIN main_bioproject mb ON b.bioproject_id = mb.bioproject_id
                WHERE mb.bioproject_name = 'CBP'
            """
            cursor.execute(query)
            result = cursor.fetchall()

        df = pd.DataFrame(result, columns=[
            "gca",
            "lowest_taxon_id",
            "scientific_name",
            "asm_name",
            "asm_level",
            "gb_status",
            "genebuilder"
        ])
        df["gb_status"] = df["gb_status"].fillna("not_started")
        df = df[~df["asm_name"].str.contains("alternate", case=False, na=False)]
        df = df.drop(columns=["asm_name"])
        df = df[~df["asm_level"].str.lower().isin(["contig", "scaffold"])]
        return df.to_dict(orient="records")

    except Exception as e:
        logging.error(f"Error fetching annotation counts: {e}")
        return []

def get_hprc():
    """Returns GCA info"""
    try:
        with get_db_connection("meta") as conn:
            cursor = conn.cursor()
            query = """
                SELECT CONCAT(a.gca_chain, '.', a.gca_version) AS gca, a.lowest_taxon_id, s.scientific_name, a.asm_name,
                        a.asm_level, gb.gb_status, gb.genebuilder
                FROM assembly a
                LEFT JOIN genebuild gb on gb.assembly_id = a.assembly_id
                LEFT JOIN bioproject b on a.assembly_id = b.assembly_id
                LEFT JOIN species s ON a.lowest_taxon_id = s.lowest_taxon_id
                LEFT JOIN group_assembly g ON a.assembly_id = g.assembly_id
                JOIN main_bioproject mb ON b.bioproject_id = mb.bioproject_id
                WHERE mb.bioproject_name = 'HPRC'
            """
            cursor.execute(query)
            result = cursor.fetchall()

        df = pd.DataFrame(result, columns=[
            "gca",
            "lowest_taxon_id",
            "scientific_name",
            "asm_name",
            "asm_level",
            "gb_status",
            "genebuilder"
        ])
        df["gb_status"] = df["gb_status"].fillna("not_started")
        df = df[~df["asm_name"].str.contains("alternate", case=False, na=False)]
        df = df.drop(columns=["asm_name"])
        df = df[~df["asm_level"].str.lower().isin(["contig", "scaffold"])]

        return df.to_dict(orient="records")

    except Exception as e:
        logging.error(f"Error fetching annotation counts: {e}")
        return []