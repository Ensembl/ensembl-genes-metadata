import json
import logging
import pandas as pd
import datetime
from fastapi import HTTPException
from metadata_app.backend.app.core.database import get_db_connection
from metadata_app.backend.app.services.taxonomy_service import get_descendant_taxa

def load_bioproject_mapping():
    """Hardcoded path for clade settings."""
    json_file = "data/bioproject_mapping.json"
    with open(json_file, "r") as f:
        logging.info("Loading bioproject mapping json file.")
        return json.load(f)


def query_meta_registry(annotation_date, taxon_id, bioproject_id, release_type):
    """Checks if each annotated assembly is the latest available version."""
    try:
        # Connect to database
        with get_db_connection("meta") as conn:
            cursor = conn.cursor()

            # Build dynamic SQL filtering
            conditions = []
            parameters = []

            if bioproject_id:
                cursor.execute("SELECT bioproject_name FROM main_bioproject")
                known_names = {row["bioproject_name"] for row in cursor.fetchall()}

                bioproject_name = [bp for bp in bioproject_id if bp in known_names]
                bioproject_ids = [bp for bp in bioproject_id if bp not in known_names]

                if bioproject_name:
                    conditions.append(f"mb.bioproject_name IN ({','.join(['%s'] * len(bioproject_name))})")
                    parameters.extend(bioproject_name)
                    logging.info(f"Filtering by BioProject names: {', '.join(bioproject_name)}")
                if bioproject_ids:
                    conditions.append(f"b.bioproject_id IN ({','.join(['%s'] * len(bioproject_ids))})")
                    parameters.extend(bioproject_ids)
                    logging.info(f"Filtering by BioProject IDs: {', '.join(bioproject_ids)}")

            if release_type:
                conditions.append(f"gb.release_type IN ({','.join(['%s'] * len(release_type))})")
                parameters.extend(release_type)
                logging.info(f"Filtering by Release Type: {', '.join(release_type)}")

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
                        detail=f"No descendant taxa found for any of the provided Taxon IDs: {', '.join(map(str, taxon_id))}"
                    )

                conditions.append(f"s.lowest_taxon_id IN ({','.join(['%s'] * len(all_descendant_taxa))})")
                parameters.extend(all_descendant_taxa)
                logging.info(f"Filtering by lowest taxon IDs: {', '.join(str(id) for id in all_descendant_taxa)}")

            if annotation_date:
                logging.info(f"Retrieving annotation for annotation date {annotation_date}.")
                if isinstance(annotation_date, pd.Timestamp):
                    annotation_date = annotation_date.strftime('%Y-%m-%d')
                elif isinstance(annotation_date, (datetime.date, datetime.datetime)):
                    annotation_date = annotation_date.strftime('%Y-%m-%d')
                conditions.append("gb.date_completed_beta >= %s")
                parameters.append(annotation_date)

            # If there are conditions, join them with AND; otherwise, select all
            where_clause = " WHERE " + " AND ".join(conditions) if conditions else ""

            meta_query = f"""
                SELECT b.bioproject_id, mb.bioproject_name AS associated_project, g.group_name, CONCAT(a.gca_chain, '.', a.gca_version) AS gca, a.lowest_taxon_id, 
                        gb.gb_status, gb.genebuilder, gb.annotation_source, gb.annotation_method, gb.release_type, 
                        gb.date_completed_beta, gb.release_date, gb.release_date_beta, gb.release_version_beta, gb.release_version,
                        s.scientific_name, s.common_name
                FROM genebuild gb
                LEFT JOIN assembly a on gb.assembly_id = a.assembly_id
                LEFT JOIN bioproject b on a.assembly_id = b.assembly_id
                LEFT JOIN species s ON a.lowest_taxon_id = s.lowest_taxon_id
                LEFT JOIN group_assembly g ON a.assembly_id = g.assembly_id
                LEFT JOIN main_bioproject mb ON b.bioproject_id = mb.bioproject_id
                {where_clause};                
            """

            cursor.execute(meta_query, parameters)
            results = cursor.fetchall()
            if not results:
                raise HTTPException(
                    status_code=404,
                    detail="No annotations found matching the specified criteria."
                )
            logging.info(f"Query returned {len(results)} annotations.")

        df_meta_genebuild = pd.DataFrame(results)
        df_meta_genebuild.drop_duplicates(subset=["gca"], inplace=True)
        logging.info(f"Retrieved records metadata table: {df_meta_genebuild.shape}")

        return df_meta_genebuild


    except Exception as e:
        logging.error(f"Unexpected error in query_meta_registry: {e}", exc_info=True)
        raise HTTPException(
            status_code=500,
            detail=f"Internal server error occurred while processing annotations: {str(e)}"
        )

def get_annotation_info_beta(df_meta_genebuild):
    """
        Fetch all annotations and their metrics, filter results based on given criteria.
        Args:
            df_meta_genebuild: Annotations table

        Returns:
            anno_wide: DataFrame with additional information
        """
    # Extract GCA accession before the version suffix

    gca_list = df_meta_genebuild['gca'].unique().tolist()
    placeholders = ', '.join(['%s'] * len(gca_list))
    logging.info(f"Retrieving annotations for GCAs: {', '.join(gca_list)}")



    try:
        # Connect to database
        with get_db_connection("prod") as conn:
            cursor = conn.cursor()

            logging.info(f"Connected to database.")

            query = f"""
                    SELECT da.value, da.attribute_id, a.accession AS gca, r.release_date AS ensembl_release_date
                    FROM dataset_attribute da
                    JOIN dataset d ON da.dataset_id = d.dataset_id
                    JOIN genome_dataset gd ON d.dataset_id = gd.dataset_id
                    JOIN genome g ON gd.genome_id = g.genome_id
                    JOIN assembly a ON g.assembly_id = a.assembly_id
                    LEFT JOIN ensembl_release r ON r.release_id = gd.release_id
                    WHERE da.attribute_id IN (34, 37, 25, 183, 31, 40, 42, 44, 48, 170, 56, 212)
                    AND a.accession IN ({placeholders})
                """

            cursor.execute(query, gca_list)
            results = cursor.fetchall()
            # Check if we have any results from the main query
            if not results:
                raise HTTPException(
                    status_code=404,
                    detail="No annotations found in production DB."
                )
            logging.info(f"Query returned {len(results)} records.")

        df_prod = pd.DataFrame(results)
        logging.info(f"Total records fetched from production db: {len(df_prod)}")

        # Directly rename attribute_id column values
        df_prod['attribute_id'] = df_prod['attribute_id'].replace({
            34: "last_geneset_update",
            37: "genebuild_method",
            212: "busco_protein",
            25: "average_exon_length",
            183: "average_sequence_length",
            31: "coding_transcripts_per_gene",
            40: "nc_average_exon_length",
            42: "nc_average_sequence_length",
            44: "nc_long_non_coding_genes",
            48: "nc_small_non_coding_genes",
            170: "nc_total_exons",
            56: "ps_average_sequence_length",
        })

        # Fill missing 'ensembl_release_date' with a placeholder value
        df_prod['ensembl_release_date'] = df_prod['ensembl_release_date'].fillna('Not yet released')

        # Pivot to wide format without losing records
        df_pivoted = df_prod.pivot_table(
            index=["ensembl_release_date", "gca"],
            columns="attribute_id",
            values="value",
            aggfunc="first"
        ).reset_index()

        logging.info(f"Pivoted DataFrame shape before filtering: {df_pivoted.shape}")

        anno_wide = pd.merge(df_meta_genebuild, df_pivoted, on='gca', how='left')

        return anno_wide



    except HTTPException:

        # Re-raise HTTPExceptions as they are already properly formatted

        raise

    except Exception as e:
        logging.error(f"Unexpected error in get_annotation_info_beta: {e}", exc_info=True)
        raise HTTPException(
            status_code=500,
            detail=f"Internal server error occurred while processing assemblies: {str(e)}"
        )


def check_if_gca_is_latest_annotated(anno_wide):
    taxon_id_list = anno_wide['lowest_taxon_id'].unique().tolist()
    placeholders = ', '.join(['%s'] * len(taxon_id_list))
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
        update_df = pd.DataFrame(results, columns=['full_gca', 'lowest_taxon_id'])
        update_df['version'] = update_df['full_gca'].str.extract(r'GCA_\d+\.(\d+)').astype(float)
        update_df['gca_root'] = update_df['full_gca'].str.replace(r'\.\d+$', '', regex=True)

        # Keep only the latest version for each root GCA
        latest_versions = (
            update_df.sort_values('version', ascending=False)
            .drop_duplicates('gca_root', keep='first')
            .rename(columns={'version': 'latest_version'})
            [['gca_root', 'latest_version']]
        )

        # Prepare the annotation DataFrame
        ann = anno_wide.copy()
        ann['version'] = ann['gca'].str.extract(r'GCA_\d+\.(\d+)').astype(float)
        ann['gca_root'] = ann['gca'].str.replace(r'\.\d+$', '', regex=True)

        # Merge to get the latest version info
        merged = ann.merge(latest_versions, on='gca_root', how='left')

        # Compare versions
        def check_latest_annotated(row):
            if pd.isna(row['latest_version']):
                return 'Yes, low quality assembly version'
            return 'Yes' if row['version'] == row['latest_version'] else 'No'

        merged['annotated_version'] = merged['version']
        merged['assembly_version'] = merged['latest_version']
        merged['latest_annotated'] = merged.apply(check_latest_annotated, axis=1)

        return merged


    except Exception as e:
        logging.error(f"Unexpected error in check_if_gca_is_latest_annotated: {e}", exc_info=True)
        raise HTTPException(
            status_code=500,
            detail=f"Internal server error occurred while processing annotations: {str(e)}"

        )


def generate_tables(annotation_date, taxon_id, bioproject_id, release_type):
    logging.info(f"Generating tables for annotation date: {annotation_date}, taxon_id: {taxon_id}, bioproject_id: {bioproject_id}, release_type: {release_type}")
    try:
        df_meta_genebuild= query_meta_registry(annotation_date, taxon_id, bioproject_id, release_type)
    except HTTPException:
        logging.error("HTTPException raised during annotation filtering")
        raise
    except Exception as e:
        logging.error("Unexpected error occurred during annotations filtering", exc_info=True)
        raise HTTPException(status_code=500, detail="An unexpected error occurred.")

    logging.info(f"Adding additional info from beta prod server")
    logging.info(f"Original df_meta_genebuild: {df_meta_genebuild.shape}")

    # Split the dataframe into beta and non-beta
    df_meta_genebuild_beta = df_meta_genebuild[df_meta_genebuild["release_type"] == "beta"].copy()
    logging.info(f"Beta subset: {df_meta_genebuild_beta.shape}")

    df_meta_genebuild_other = df_meta_genebuild[df_meta_genebuild["release_type"] != "beta"].copy()
    logging.info(f"Other subset: {df_meta_genebuild_other.shape}")

    if not df_meta_genebuild_beta.empty:
        df_meta_genebuild_beta_updated = get_annotation_info_beta(df_meta_genebuild_beta)

        # If the function fails and returns None, fallback to original beta with "Error"
        if df_meta_genebuild_beta_updated is None:
            df_meta_genebuild_beta["latest_annotated"] = "Error"
            df_meta_genebuild_beta_updated = df_meta_genebuild_beta
    else:
        df_meta_genebuild_beta_updated = pd.DataFrame(columns=df_meta_genebuild.columns.tolist() + ["latest_annotated"])
    logging.info(f"Updated beta: {df_meta_genebuild_beta_updated.shape}")

    # Concatenate the updated beta subset with the rest
    anno_wide_pre_check = pd.concat([df_meta_genebuild_beta_updated, df_meta_genebuild_other], ignore_index=True)
    logging.info(f"Final combined: {anno_wide_pre_check.shape}")
    logging.info(f"After beta info check: {anno_wide_pre_check.shape}")


    logging.info(f"Checking if annotation is the latest GCA version")
    anno_wide = check_if_gca_is_latest_annotated(anno_wide_pre_check)
    logging.info(f"After latest annotated check: {anno_wide.shape}")

    # Create the FTP URL using the scientific_name, replacing spaces with underscores
    logging.info("Generating FTP paths.")
    anno_wide['ftp'] = anno_wide.apply(
        lambda
            row: f"https://ftp.ebi.ac.uk/pub/ensemblorganisms/{row['scientific_name'].replace(' ', '_')}/{row['gca']}/"
        if pd.notnull(row['scientific_name']) and pd.notnull(row['gca']) and pd.notnull(
            row.get('ensembl_release_date')) else None,
        axis=1)


    #filtered_df = filtered_df.drop(columns=['year', 'gca', 'version'])
    #df_info_result = df_info_result.drop(columns=['year', 'version', 'gca_latest'])
    anno_wide = anno_wide.drop_duplicates(subset='gca', keep='first')
    # Create main display table
    anno_main = anno_wide[
        ['bioproject_id', 'associated_project', 'gca', 'scientific_name', 'date_completed_beta',
         'release_date_beta', 'lowest_taxon_id', 'gb_status', 'release_type', 'latest_annotated']
    ]

    # Transforming out of range float values that are not JSON compliant: nan
    logging.info(f"Transfroming Out of range float values that are not JSON compliant")
    anno_main = anno_main.apply(lambda col: col.fillna("") if col.dtype == "object" else col)
    anno_wide = anno_wide.apply(lambda col: col.fillna("") if col.dtype == "object" else col)

    return anno_wide, anno_main