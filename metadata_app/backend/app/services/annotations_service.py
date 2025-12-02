import logging
import numpy as np
import pandas as pd
import datetime
from fastapi import HTTPException
from metadata_app.backend.app.core.database import get_db_connection
from metadata_app.backend.app.services.taxonomy_service import get_descendant_taxa, load_clade_data, \
    assign_clade_and_species


def query_meta_registry(annotation_date, taxon_id, bioproject_id, group_name):
    """Checks if each annotated assembly is the latest available version."""
    try:
        # Connect to database
        with get_db_connection("meta") as conn:
            cursor = conn.cursor()

            # Validate BioProject IDs if provided
            if bioproject_id:
                cursor.execute("SELECT DISTINCT bioproject_id FROM bioproject;")
                valid_bioprojects = {row['bioproject_id'] for row in cursor.fetchall()}
                invalid_bioprojects = set(bioproject_id) - valid_bioprojects
                if invalid_bioprojects:
                    raise HTTPException(
                        status_code=400,
                        detail=f"The following BioProject IDs were not found: {', '.join(invalid_bioprojects)}"
                    )

            # Build dynamic SQL filtering
            conditions = []
            parameters = []

            if bioproject_id:
                conditions.append(f"b.bioproject_id IN ({','.join(['%s'] * len(bioproject_id))})")
                parameters.extend(bioproject_id)
                logging.info(f"Filtering by BioProject IDs: {', '.join(bioproject_id)}")

            if group_name:
                conditions.append("g.group_name = %s")
                parameters.append(group_name)
                logging.info(f"Filtering by group name: {group_name}")

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
                conditions.append("gb.date_status_update >= %s")
                parameters.append(annotation_date)

            # If there are conditions, join them with AND; otherwise, select all
            where_clause = " WHERE " + " AND ".join(conditions) if conditions else ""

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
                    am.protein_busco,
                    am.protein_busco_lineage,
                    am.protein_busco_version
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
                             MAX(CASE WHEN metrics_name='genebuild_busco' THEN metrics_value END) AS protein_busco,
                             MAX(CASE WHEN metrics_name='genebuild_busco_dataset' THEN metrics_value END) AS protein_busco_lineage,
                             MAX(CASE WHEN metrics_name='genebuild_busco_version' THEN metrics_value END) AS protein_busco_version
                      FROM annotation_metrics
                      GROUP BY genebuild_status_id
                ) am ON gb.genebuild_status_id = am.genebuild_status_id
                {where_clause}
                AND gb.last_attempt = 1
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
                    detail="No annotations found matching the specified criteria."
                )
            logging.info(f"Query returned {len(results)} lines.")

            # Get taxonomy data
            lowest_taxon_ids = {row['lowest_taxon_id'] for row in results if
                                'lowest_taxon_id' in row and row['lowest_taxon_id'] is not None}
            logging.debug(f"Collected lowest taxon IDs {print(lowest_taxon_ids)}")

            if not lowest_taxon_ids:
                # No results or no taxon IDs found
                raise HTTPException(status_code=404, detail="No valid taxon IDs found in the results.")

            # Fetch all taxonomy data for the collected lowest_taxon_ids
            taxonomy_query = """
                            SELECT lowest_taxon_id, taxon_class_id, taxon_class
                            FROM taxonomy
                            WHERE lowest_taxon_id IN ({})
                            ORDER BY FIELD(taxon_class, 'species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom');
                        """.format(','.join(['%s'] * len(lowest_taxon_ids)))

            cursor.execute(taxonomy_query, tuple(lowest_taxon_ids))
            taxonomy_results = cursor.fetchall()
            logging.info(f"Taxonomy Query executed successfully, retrieved {len(taxonomy_results)} results.")

            # Process taxonomy results
            taxonomy_dict = {}
            for row in taxonomy_results:
                lowest_taxon_id = row['lowest_taxon_id']
                if lowest_taxon_id not in taxonomy_dict:
                    taxonomy_dict[lowest_taxon_id] = []
                taxonomy_dict[lowest_taxon_id].append({
                    'taxon_class_id': row['taxon_class_id'],
                    'taxon_class': row['taxon_class']
                })



        df_meta_genebuild = pd.DataFrame(results)

        # Add clade, species, and genus information
        clade_data = load_clade_data()

        df_meta_genebuild[['internal_clade', 'species_taxon_id', 'genus_taxon_id']] = df_meta_genebuild[
            'lowest_taxon_id'].apply(
            lambda x: pd.Series(assign_clade_and_species(x, clade_data, taxonomy_dict))
        )

        logging.info(f"Added clade data")
        logging.info(f"Changing genus id format")
        df_meta_genebuild['genus_taxon_id'] = (
            pd.to_numeric(df_meta_genebuild['genus_taxon_id'].replace('', pd.NA), errors='coerce')
            .astype('Int64')
        )
        logging.info(f"Changed genus id format")

        logging.info(f"Retrieved records from genebuild_status table: {df_meta_genebuild.shape}")
        print(df_meta_genebuild)
        return df_meta_genebuild


    except Exception as e:
        logging.error(f"Unexpected error in query_meta_registry: {e}", exc_info=True)
        raise HTTPException(
            status_code=500,
            detail=f"Internal server error occurred while processing annotations: {str(e)}"
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


def generate_tables(annotation_date, taxon_id, bioproject_id, group_name):
    logging.info(f"Generating tables for annotation date: {annotation_date}, taxon_id: {taxon_id}, bioproject_id: {bioproject_id}")
    try:
        df_meta_genebuild= query_meta_registry(annotation_date, taxon_id, bioproject_id, group_name)
    except HTTPException:
        logging.error("HTTPException raised during annotation filtering")
        raise
    except Exception as e:
        logging.error("Unexpected error occurred during annotations filtering", exc_info=True)
        raise HTTPException(status_code=500, detail="An unexpected error occurred.")



    logging.info(f"Checking if annotation is the latest GCA version")
    anno_wide = check_if_gca_is_latest_annotated(df_meta_genebuild)
    logging.info(f"After latest annotated check: {anno_wide.shape}")

    # Create the FTP URL using the scientific_name, replacing spaces with underscores
    logging.info("Generating FTP paths.")
    anno_wide['ftp'] = anno_wide.apply(
        lambda
            row: f"https://ftp.ebi.ac.uk/pub/ensemblorganisms/{row['scientific_name'].replace(' ', '_')}/{row['gca']}/"
        if pd.notnull(row['scientific_name']) and pd.notnull(row['gca']) and pd.notnull(
            row.get('release_date')) else None,
        axis=1)


    #filtered_df = filtered_df.drop(columns=['year', 'gca', 'version'])
    #df_info_result = df_info_result.drop(columns=['year', 'version', 'gca_latest'])
    anno_wide = anno_wide.drop_duplicates(subset='gca', keep='first')
    # Create main display table
    anno_main = anno_wide[
        ['bioproject_id', 'associated_project', 'gca', 'scientific_name', 'last_genebuild_update',
         'release_date', 'lowest_taxon_id', 'gb_status', 'latest_annotated']
    ]

    # Transforming out of range float values that are not JSON compliant: nan
    logging.info(f"Transfroming Out of range float values that are not JSON compliant")
    anno_main = anno_main.apply(lambda col: col.fillna("") if col.dtype == "object" else col)
    anno_wide = anno_wide.apply(lambda col: col.fillna("") if col.dtype == "object" else col)

    return anno_wide, anno_main