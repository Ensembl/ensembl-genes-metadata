import pandas as pd
import logging
from metadata_app.backend.app.services.annotations_service import generate_tables
from metadata_app.backend.app.services.taxonomy_service import (
    get_descendant_taxa,
    assign_clade_and_species,
    load_clade_data,
)
from fastapi import HTTPException
from metadata_app.backend.app.core.database import get_db_connection


def generate_report(end_date, start_date, group_name, taxon_id, bioproject_id):
    anno_wide, anno_main = generate_tables(
        group_name=group_name,
        taxon_id=taxon_id,
        bioproject_id=bioproject_id,
        annotation_date=None,
    )
    # Ensure date fields are datetime
    for col in ["last_genebuild_update", "release_date", "date_status_update"]:
        if col in anno_wide.columns:
            anno_wide[col] = pd.to_datetime(anno_wide[col], errors="coerce")

    if end_date:
        end_date = pd.to_datetime(end_date)
        anno_wide = anno_wide[anno_wide["date_status_update"] <= end_date]

    if start_date:
        start_date = pd.to_datetime(start_date)
        anno_wide = anno_wide[anno_wide["date_status_update"] >= start_date]

    # Create tables for charts
    number_of_annotations_raw = (
        anno_wide[["gca", "gb_status"]]
        .groupby("gb_status")
        .size()
        .reset_index(name="count")
    )

    # Transform into desired format
    number_of_annotations = [
        {"gb_status": row["gb_status"], "count": row["count"]}
        for _, row in number_of_annotations_raw.iterrows()
    ]

    method_report = (
        anno_wide[["gca", "annotation_method"]]
        .groupby("annotation_method")
        .size()
        .reset_index(name="count")
    )

    num_unique_taxa = anno_wide["lowest_taxon_id"].nunique()
    top_3_taxa = (
        anno_wide.groupby(["scientific_name"])
        .size()
        .reset_index(name="count")
        .sort_values(by="count", ascending=False)
        .head(3)
    )

    project_report = (
        anno_wide[["gca", "associated_project"]]
        .groupby("associated_project")
        .size()
        .reset_index(name="count")
    )

    if "protein_busco" in anno_wide.columns:
        extracted = anno_wide["protein_busco"].str.extract(r"C:(\d+\.\d+)%")[0]

        # Convert to float, but invalid values → None instead of NaN
        extracted = pd.to_numeric(extracted, errors="coerce")
        # Store cleaned series back in the DataFrame (optional)
        anno_wide["busco_complete"] = extracted
        # Compute average only on valid numbers
        valid = extracted.dropna()

        average_busco = valid.mean() if not valid.empty else "Not available"
    else:
        anno_wide["protein_busco"] = "Not available"
        average_busco = "Not available"

    main_report = anno_wide[
        [
            "associated_project",
            "gca",
            "genebuilder",
            "gb_status",
            "ftp",
            "latest_annotated",
            "protein_busco",
            "last_genebuild_update",
            "release_date",
        ]
    ]
    logging.info(f"BUSCo {average_busco}")

    # Get clade information
    # Add clade, species, and genus information

    # Get taxonomy data
    lowest_taxon_ids = set(anno_wide["lowest_taxon_id"].dropna().unique())
    logging.debug(f"Collected lowest taxon IDs: {lowest_taxon_ids}")

    if not lowest_taxon_ids:
        # No results or no taxon IDs found
        raise HTTPException(
            status_code=404, detail="No valid taxon IDs found in the results."
        )

    # Fetch all taxonomy data for the collected lowest_taxon_ids
    try:
        with get_db_connection("meta") as conn:
            cursor = conn.cursor()
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
    except Exception as e:
        logging.error(f"Taxonomy query failed: {e}")
        raise HTTPException(status_code=500, detail="Error fetching taxonomy data")

    # Process taxonomy results
    taxonomy_dict = {}
    for row in taxonomy_results:
        lowest_taxon_id = row["lowest_taxon_id"]
        if lowest_taxon_id not in taxonomy_dict:
            taxonomy_dict[lowest_taxon_id] = []
        taxonomy_dict[lowest_taxon_id].append(
            {"taxon_class_id": row["taxon_class_id"], "taxon_class": row["taxon_class"]}
        )
    clade_data = load_clade_data()

    anno_wide[["internal_clade", "species_taxon_id", "genus_taxon_id", "pipeline"]] = anno_wide[
        "lowest_taxon_id"
    ].apply(lambda x: pd.Series(assign_clade_and_species(x, clade_data, taxonomy_dict)))

    logging.info("Added clade data")
    logging.info("Changing genus id format")
    anno_wide["genus_taxon_id"] = pd.to_numeric(
        anno_wide["genus_taxon_id"].replace("", pd.NA), errors="coerce"
    ).astype("Int64")
    logging.info("Changed genus id format")

    # Create clade summary
    clade_group = (
        anno_wide[["gca", "internal_clade"]]
        .groupby("internal_clade")
        .size()
        .reset_index(name="count")
    )
    logging.info("Created clade group summary")

    # Create clade summary for live annotations
    live_anno = anno_wide[anno_wide["gb_status"] == "live"]
    clade_group_live = (
        live_anno[["gca", "internal_clade"]]
        .groupby("internal_clade")
        .size()
        .reset_index(name="count")
    )
    logging.info("Created live clade group summary")

    # Transforming out of range float values that are not JSON compliant: nan
    logging.info(f"Transfroming Out of range float values that are not JSON compliant")
    anno_wide = anno_wide.apply(
        lambda col: col.fillna("") if col.dtype == "object" else col
    )
    method_report = method_report.apply(
        lambda col: col.fillna("") if col.dtype == "object" else col
    )
    top_3_taxa = top_3_taxa.apply(
        lambda col: col.fillna("") if col.dtype == "object" else col
    )
    project_report = project_report.apply(
        lambda col: col.fillna("") if col.dtype == "object" else col
    )
    main_report = main_report.apply(
        lambda col: col.fillna("") if col.dtype == "object" else col
    )
    clade_group = clade_group.apply(
        lambda col: col.fillna("") if col.dtype == "object" else col
    )
    clade_group_live = clade_group_live.apply(
        lambda col: col.fillna("") if col.dtype == "object" else col
    )

    return (
        anno_wide,
        number_of_annotations,
        method_report,
        num_unique_taxa,
        top_3_taxa,
        project_report,
        average_busco,
        main_report,
        clade_group,
        clade_group_live,
    )
