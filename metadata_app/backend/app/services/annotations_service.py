"""
This module provides services for handling annotations-related operations.
"""

import logging
import numpy as np
import pandas as pd
import datetime
import re
from fastapi import HTTPException
from metadata_app.backend.app.core.database import get_db_connection
from metadata_app.backend.app.services.taxonomy_service import (
    get_descendant_taxa,
    load_clade_data,
    assign_clade_and_species,
)

ENSEMBL_ORGANISMS_FTP_BASE_URL = "https://ftp.ebi.ac.uk/pub/ensemblorganisms"


def _format_assembly_accession_path(gca_accession):
    """
    Formats the given GCA accession into a path structure.

    Args:
        gca_accession: The GCA accession to be formatted.

    Returns:
        str: The formatted path or None if the input is invalid.
    """
    accession = str(gca_accession).strip()
    match = re.fullmatch(
        r"(?P<prefix>GC[AF])_(?P<digits>\d+)\.(?P<version>\d+)", accession
    )
    if not match:
        return None

    digits = match.group("digits")
    if len(digits) % 3 != 0:
        return None

    grouped_digits = "/".join(
        digits[index : index + 3] for index in range(0, len(digits), 3)
    )
    return f"{match.group('prefix')}/{grouped_digits}/{match.group('version')}"


def _format_annotation_date_path(date_value):
    if pd.isna(date_value):
        return None

    timestamp = pd.to_datetime(date_value, errors="coerce")
    if pd.isna(timestamp):
        return None

    return timestamp.strftime("%Y_%m")


def _format_annotation_provider_path(annotation_source):
    """
    Formats the annotation source to determine the appropriate provider path.

    This function processes the given annotation source and standardizes it
    to one of the predefined annotation provider paths. The logic accounts for
    empty, null, or whitespace strings, converting them to the "ensembl" path.
    A valid, non-empty source other than "ensembl" will default to "community".

    Args:
        annotation_source: The input source of the annotation. May be a string,
            null, or NaN value.

    Returns:
        str: A string representing the formatted annotation provider path. It
        will return either "ensembl" or "community" based on the input source.
    """
    if pd.isna(annotation_source):
        return "ensembl"

    source = str(annotation_source).strip().lower()
    if not source:
        return "ensembl"

    return "ensembl" if source == "ensembl" else "community"


def build_ensemblorganisms_ftp_url(row):
    """
    Builds the FTP URL for Ensembl organisms based on the provided row data.
    """
    accession_path = _format_assembly_accession_path(row.get("gca"))
    annotation_date = _format_annotation_date_path(row.get("last_genebuild_update"))
    provider = _format_annotation_provider_path(row.get("annotation_source"))

    if not accession_path or not annotation_date or pd.isna(row.get("release_date")):
        return None

    return (
        f"{ENSEMBL_ORGANISMS_FTP_BASE_URL}/"
        f"{accession_path}/{provider}/{annotation_date}/"
    )


def query_meta_registry(annotation_date, taxon_id, bioproject_id, group_name, gca):
    """Get annotataions"""
    try:
        # Connect to database
        with get_db_connection("meta") as conn:
            cursor = conn.cursor()

            if isinstance(bioproject_id, str):
                bioproject_id = [bioproject_id]
            if isinstance(group_name, str):
                group_name = [group_name]

            # Validate BioProject IDs if provided
            if bioproject_id:
                cursor.execute("SELECT DISTINCT bioproject_id FROM bioproject;")
                valid_bioprojects = {row["bioproject_id"] for row in cursor.fetchall()}
                invalid_bioprojects = set(bioproject_id) - valid_bioprojects
                if invalid_bioprojects:
                    raise HTTPException(
                        status_code=400,
                        detail=f"The following BioProject IDs were not found: {', '.join(invalid_bioprojects)}",
                    )

            # Build dynamic SQL filtering
            conditions = []
            parameters = []
            project_conditions = []

            if bioproject_id:
                project_conditions.append(
                    f"b.bioproject_id IN ({','.join(['%s'] * len(bioproject_id))})"
                )
                parameters.extend(bioproject_id)
                logging.info(f"Filtering by BioProject IDs: {', '.join(bioproject_id)}")

            if group_name:
                project_conditions.append(
                    f"g.group_name IN ({','.join(['%s'] * len(group_name))})"
                )
                parameters.extend(group_name)
                logging.info(f"Filtering by group name: {', '.join(group_name)}")

            if len(project_conditions) == 1:
                conditions.append(project_conditions[0])
            elif len(project_conditions) > 1:
                conditions.append(f"({' OR '.join(project_conditions)})")

            if taxon_id:
                all_descendant_taxa = set()
                for tax_id in taxon_id:
                    descendant_taxa = get_descendant_taxa(tax_id)
                    if not descendant_taxa:
                        logging.warning(f"No descendants found for taxon ID {tax_id}")
                    all_descendant_taxa.update(descendant_taxa)

                if not all_descendant_taxa:
                    raise HTTPException(
                        status_code=400,
                        detail=f"No descendant taxa found for any of the provided Taxon IDs: {', '.join(map(str, taxon_id))}",
                    )

                conditions.append(
                    f"s.lowest_taxon_id IN ({','.join(['%s'] * len(all_descendant_taxa))})"
                )
                parameters.extend(all_descendant_taxa)
                logging.info(
                    f"Filtering by lowest taxon IDs: {', '.join(str(id) for id in all_descendant_taxa)}"
                )

            if annotation_date:
                logging.info(
                    f"Retrieving annotation for annotation date {annotation_date}."
                )
                if isinstance(annotation_date, pd.Timestamp):
                    annotation_date = annotation_date.strftime("%Y-%m-%d")
                elif isinstance(annotation_date, (datetime.date, datetime.datetime)):
                    annotation_date = annotation_date.strftime("%Y-%m-%d")
                conditions.append("gb.date_status_update >= %s")
                parameters.append(annotation_date)

            if gca:
                if isinstance(gca, str):
                    gca = [gca]
                gca_list_filter = ",".join(["%s"] * len(gca))
                conditions.append(f"gb.gca_accession IN ({gca_list_filter})")
                parameters.extend(gca)
                logging.info(f"Filtering by GCA: {', '.join(gca)}")

            conditions.append("gb.last_attempt = 1")
            where_clause = " WHERE " + " AND ".join(conditions)

            meta_query = f"""
                SELECT 
                    b.bioproject_id,
                    mb.bioproject_name AS associated_project,
                    g.group_name,
                    CONCAT(a.gca_chain, '.', a.gca_version) AS gca,
                    a.lowest_taxon_id,
                    gb.gb_status,
                    gb.genebuilder,
                    gb.annotation_source,
                    gb.annotation_method,
                    gb.date_started,
                    gb.release_date,
                    gb.date_status_update,
                    gb.last_genebuild_update,
                    s.scientific_name,
                    s.common_name,
                    asm.assembly_busco,
                    asm.assembly_busco_lineage,
                    asm.assembly_busco_version,
                    am.protein_busco,
                    am.protein_busco_lineage,
                    am.protein_busco_version,
                    am.coding_genes
                FROM genebuild_status gb
                LEFT JOIN assembly a ON gb.assembly_id = a.assembly_id
                LEFT JOIN bioproject b ON a.assembly_id = b.assembly_id
                LEFT JOIN species s ON a.lowest_taxon_id = s.lowest_taxon_id
                LEFT JOIN custom_group g
                    ON (
                         (g.group_type = 'taxon' AND a.lowest_taxon_id = g.item)
                         OR
                         (g.group_type = 'assembly' AND a.gca_chain = g.item)
                       )
                LEFT JOIN main_bioproject mb ON b.bioproject_id = mb.bioproject_id
                LEFT JOIN (
                      SELECT genebuild_status_id,
                             MAX(CASE WHEN metrics_name='genebuild.busco' THEN metrics_value END) AS protein_busco,
                             MAX(CASE WHEN metrics_name='genebuild.busco_dataset' THEN metrics_value END) AS protein_busco_lineage,
                             MAX(CASE WHEN metrics_name='genebuild.busco_version' THEN metrics_value END) AS protein_busco_version,
                             MAX(CASE WHEN metrics_name='genebuild.stats.coding_genes' THEN metrics_value END) AS coding_genes
                      FROM annotation_metrics
                      GROUP BY genebuild_status_id
                ) am ON gb.genebuild_status_id = am.genebuild_status_id
                LEFT JOIN (
                      SELECT assembly_id,
                             MAX(CASE WHEN metrics_name='assembly.busco' THEN metrics_value END) AS assembly_busco,
                             MAX(CASE WHEN metrics_name='assembly.busco_dataset' THEN metrics_value END) AS assembly_busco_lineage,
                             MAX(CASE WHEN metrics_name='assembly.busco_version' THEN metrics_value END) AS assembly_busco_version
                      FROM assembly_metrics
                      GROUP BY assembly_id
                ) asm ON a.assembly_id = asm.assembly_id
                {where_clause}
                GROUP BY
                    b.bioproject_id,
                    mb.bioproject_name,
                    g.group_name,
                    a.gca_chain,
                    a.gca_version,
                    a.lowest_taxon_id,
                    gb.gb_status,
                    gb.genebuilder,
                    gb.annotation_source,
                    gb.annotation_method,
                    gb.date_started,
                    gb.release_date,
                    gb.date_status_update,
                    gb.last_genebuild_update,
                    s.scientific_name,
                    s.common_name;                
            """

            cursor.execute(meta_query, parameters)
            results = cursor.fetchall()
            if not results:
                raise HTTPException(
                    status_code=404,
                    detail="No annotations found matching the specified criteria.",
                )
            logging.info(f"Query returned {len(results)} lines.")

            # Get taxonomy data
            lowest_taxon_ids = {
                row["lowest_taxon_id"]
                for row in results
                if "lowest_taxon_id" in row and row["lowest_taxon_id"] is not None
            }
            logging.debug(f"Collected lowest taxon IDs {print(lowest_taxon_ids)}")

            if not lowest_taxon_ids:
                # No results or no taxon IDs found
                raise HTTPException(
                    status_code=404, detail="No valid taxon IDs found in the results."
                )

            # Fetch all taxonomy data for the collected lowest_taxon_ids
            taxonomy_query = """
                            SELECT lowest_taxon_id, taxon_class_id, taxon_class
                            FROM taxonomy
                            WHERE lowest_taxon_id IN ({})
                            ORDER BY FIELD(taxon_class, 'species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom');
                        """.format(",".join(["%s"] * len(lowest_taxon_ids)))

            cursor.execute(taxonomy_query, tuple(lowest_taxon_ids))
            taxonomy_results = cursor.fetchall()
            logging.info(
                f"Taxonomy Query executed successfully, retrieved {len(taxonomy_results)} results."
            )

            # Process taxonomy results
            taxonomy_dict = {}
            for row in taxonomy_results:
                lowest_taxon_id = row["lowest_taxon_id"]
                if lowest_taxon_id not in taxonomy_dict:
                    taxonomy_dict[lowest_taxon_id] = []
                taxonomy_dict[lowest_taxon_id].append(
                    {
                        "taxon_class_id": row["taxon_class_id"],
                        "taxon_class": row["taxon_class"],
                    }
                )

        df_meta_genebuild = pd.DataFrame(results)

        # Add clade, species, and genus information
        clade_data = load_clade_data()

        df_meta_genebuild[
            ["internal_clade", "species_taxon_id", "genus_taxon_id", "pipeline"]
        ] = df_meta_genebuild["lowest_taxon_id"].apply(
            lambda x: pd.Series(assign_clade_and_species(x, clade_data, taxonomy_dict))
        )

        logging.info(f"Added clade data")
        logging.info(f"Changing genus id format")
        df_meta_genebuild["genus_taxon_id"] = pd.to_numeric(
            df_meta_genebuild["genus_taxon_id"].replace("", pd.NA), errors="coerce"
        ).astype("Int64")
        logging.info(f"Changed genus id format")

        logging.info(
            f"Retrieved records from genebuild_status table: {df_meta_genebuild.shape}"
        )
        print(df_meta_genebuild)
        return df_meta_genebuild

    except Exception as e:
        logging.error(f"Unexpected error in query_meta_registry: {e}", exc_info=True)
        raise HTTPException(
            status_code=500,
            detail=f"Internal server error occurred while processing annotations: {str(e)}",
        )


def check_if_gca_is_latest_annotated(anno_wide):
    taxon_id_list = anno_wide["lowest_taxon_id"].unique().tolist()
    placeholders = ", ".join(["%s"] * len(taxon_id_list))
    logging.info(f"taxon_id list: {taxon_id_list}")

    try:
        # Connect to database
        with get_db_connection("meta") as conn:
            cursor = conn.cursor()
            update_query = f"""
                SELECT CONCAT(a.gca_chain, '.', a.gca_version) AS full_gca, a.lowest_taxon_id
                FROM assembly a
                WHERE a.lowest_taxon_id IN ({placeholders});
            """
            cursor.execute(update_query, taxon_id_list)
            results = cursor.fetchall()

        # Convert to DataFrame
        update_df = pd.DataFrame(results, columns=["full_gca", "lowest_taxon_id"])
        update_df["version"] = (
            update_df["full_gca"].str.extract(r"GCA_\d+\.(\d+)").astype(float)
        )
        update_df["gca_root"] = update_df["full_gca"].str.replace(
            r"\.\d+$", "", regex=True
        )

        # Keep only the latest version for each root GCA
        latest_versions = (
            update_df.sort_values("version", ascending=False)
            .drop_duplicates("gca_root", keep="first")
            .rename(columns={"version": "latest_version"})[
                ["gca_root", "latest_version"]
            ]
        )

        # Prepare the annotation DataFrame
        ann = anno_wide.copy()
        ann["version"] = ann["gca"].str.extract(r"GCA_\d+\.(\d+)").astype(float)
        ann["gca_root"] = ann["gca"].str.replace(r"\.\d+$", "", regex=True)

        # Merge to get the latest version info
        merged = ann.merge(latest_versions, on="gca_root", how="left")

        # Compare versions
        def check_latest_annotated(row):
            if pd.isna(row["latest_version"]):
                return "Yes, low quality assembly version"
            return "Yes" if row["version"] == row["latest_version"] else "No"

        merged["annotated_version"] = merged["version"]
        merged["assembly_version"] = merged["latest_version"]
        merged["latest_annotated"] = merged.apply(check_latest_annotated, axis=1)

        return merged

    except Exception as e:
        logging.error(
            f"Unexpected error in check_if_gca_is_latest_annotated: {e}", exc_info=True
        )
        raise HTTPException(
            status_code=500,
            detail=f"Internal server error occurred while processing annotations: {str(e)}",
        )


def generate_tables(annotation_date, taxon_id, bioproject_id, group_name, gca=None):
    logging.info(
        f"Generating tables for annotation date: {annotation_date}, taxon_id: {taxon_id}, bioproject_id: {bioproject_id}, gca: {gca}"
    )
    try:
        df_meta_genebuild = query_meta_registry(
            annotation_date, taxon_id, bioproject_id, group_name, gca
        )
    except HTTPException:
        logging.error("HTTPException raised during annotation filtering")
        raise
    except Exception as e:
        logging.error(
            "Unexpected error occurred during annotations filtering", exc_info=True
        )
        raise HTTPException(status_code=500, detail="An unexpected error occurred.")

    logging.info(f"Checking if annotation is the latest GCA version")
    anno_wide = check_if_gca_is_latest_annotated(df_meta_genebuild)
    logging.info(f"After latest annotated check: {anno_wide.shape}")

    # Create the FTP URL using the accession-based Ensembl Organisms structure.
    logging.info("Generating FTP paths.")
    anno_wide["ftp"] = anno_wide.apply(build_ensemblorganisms_ftp_url, axis=1)
    anno_project_memberships = anno_wide[
        ["gca", "associated_project", "gb_status", "date_status_update"]
    ].drop_duplicates(subset=["gca", "associated_project"])

    # filtered_df = filtered_df.drop(columns=['year', 'gca', 'version'])
    # df_info_result = df_info_result.drop(columns=['year', 'version', 'gca_latest'])
    anno_wide = anno_wide.drop_duplicates(subset="gca", keep="first")
    # Create main display table
    anno_main = anno_wide[
        [
            "bioproject_id",
            "associated_project",
            "gca",
            "scientific_name",
            "last_genebuild_update",
            "date_status_update",
            "release_date",
            "lowest_taxon_id",
            "gb_status",
            "latest_annotated",
            "protein_busco",
            "protein_busco_lineage",
            "assembly_busco",
            "assembly_busco_lineage",
        ]
    ]

    # Transforming out of range float values that are not JSON compliant: nan
    logging.info(f"Transfroming Out of range float values that are not JSON compliant")
    anno_main = anno_main.apply(
        lambda col: col.fillna("") if col.dtype == "object" else col
    )
    anno_wide = anno_wide.apply(
        lambda col: col.fillna("") if col.dtype == "object" else col
    )
    anno_project_memberships = anno_project_memberships.apply(
        lambda col: col.fillna("") if col.dtype == "object" else col
    )

    return anno_wide, anno_main, anno_project_memberships
