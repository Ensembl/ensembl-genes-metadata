import json
import logging
import pandas as pd
import datetime
from fastapi import HTTPException
from numpy.ma.extras import average

from metadata_app.backend.app.core.database import get_db_connection
from metadata_app.backend.app.services.taxonomy_service import get_descendant_taxa

def load_bioproject_mapping():
    """Hardcoded path for clade settings."""
    json_file = "data/bioproject_mapping.json"
    with open(json_file, "r") as f:
        logging.info("Loading bioproject mapping json file.")
        return json.load(f)


def query_meta_registry(start_date, end_date, group_name, taxon_id, bioproject_id, release_type):
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
                    return f"The following BioProject IDs were not found: {', '.join(invalid_bioprojects)}", None, None

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

            if release_type:
                conditions.append(f"gb.release_type IN ({','.join(['%s'] * len(release_type))})")
                parameters.extend(release_type)
                logging.info(f"Filtering by Release Type: {', '.join(release_type)}")

            if taxon_id:
                descendant_taxa = get_descendant_taxa(taxon_id)
                logging.info(f"Retrieving annotation for taxon ID {taxon_id}.")
                if not descendant_taxa:
                    return f"No descendant taxa found for Taxon ID {taxon_id}.", None

                conditions.append(f"a.lowest_taxon_id IN ({','.join(['%s'] * len(descendant_taxa))})")
                parameters.extend(descendant_taxa)

            if start_date:
                logging.info(f"Retrieving annotation for start date {start_date}.")
                if isinstance(start_date, pd.Timestamp):
                    start_date = start_date.strftime('%Y-%m-%d')
                elif isinstance(start_date, (datetime.date, datetime.datetime)):
                    start_date = start_date.strftime('%Y-%m-%d')
                conditions.append("gb.date_completed_beta >= %s")
                parameters.append(start_date)

            if end_date:
                logging.info(f"Retrieving annotation for end date {end_date}.")
                if isinstance(end_date, pd.Timestamp):
                    end_date = end_date.strftime('%Y-%m-%d')
                elif isinstance(end_date, (datetime.date, datetime.datetime)):
                    end_date = end_date.strftime('%Y-%m-%d')
                conditions.append("gb.date_completed_beta <= %s")
                parameters.append(end_date)

            # If there are conditions, join them with AND; otherwise, select all
            where_clause = " WHERE " + " AND ".join(conditions) if conditions else ""

            meta_query = f"""
                SELECT b.bioproject_id, g.group_name, CONCAT(a.gca_chain, '.', a.gca_version) AS gca, a.lowest_taxon_id, 
                        gb.gb_status, gb.genebuilder, gb.annotation_source, gb.annotation_method, gb.release_type, 
                        gb.date_completed_beta, gb.release_date, gb.release_date_beta, gb.release_version_beta, gb.release_version,
                        s.scientific_name, s.common_name
                FROM genebuild gb
                LEFT JOIN assembly a on gb.assembly_id = a.assembly_id
                LEFT JOIN bioproject b on a.assembly_id = b.assembly_id
                LEFT JOIN species s ON a.lowest_taxon_id = s.lowest_taxon_id
                LEFT JOIN group_assembly g ON a.assembly_id = g.assembly_id
                {where_clause};                
            """

            cursor.execute(meta_query, parameters)
            results = cursor.fetchall()

        df_meta_genebuild = pd.DataFrame(results)
        logging.info(f"Retrieved records: {df_meta_genebuild.shape}")


        # Load bioproject_mapping
        logging.info(f"Adding BioProject mapping")
        bioproject_mapping = load_bioproject_mapping()
        df_meta_genebuild["associated_project"] = df_meta_genebuild["bioproject_id"].map(bioproject_mapping)
        logging.info(f"Added BioProject mapping")

        return df_meta_genebuild

    except Exception as e:
        logging.error(f"Error in query_meta_registry: {str(e)}")
        return None

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


    except Exception as e:
            logging.error(f"Error in get_annotation_info_beta: {str(e)}")
            return None

def check_if_gca_is_latest_annotated(anno_wide):
    taxon_id_list = anno_wide['lowest_taxon_id'].unique().tolist()
    placeholders = ', '.join(['%s'] * len(taxon_id_list))
    logging.info(f"taxon_id list: {taxon_id_list}")
    try:
        # Connect to database
        with get_db_connection("meta") as conn:
            cursor = conn.cursor()

            update_query = f"""
                SELECT CONCAT(a.gca_chain, '.', a.gca_version) AS gca, a.lowest_taxon_id
                FROM assembly a
                WHERE a.lowest_taxon_id IN ({placeholders});
            """

            cursor.execute(update_query, taxon_id_list)

            results = cursor.fetchall()

        update_df = pd.DataFrame(results)

        # Extract version number from GCA identifier
        asm = update_df
        asm['version'] = asm['gca'].str.extract(r'GCA_\d+\.(\d+)').astype(float)
        asm['gca_latest'] = asm['gca']
        ann = anno_wide
        ann['version'] = ann['gca'].str.extract(r'GCA_\d+\.(\d+)').astype(float)
        ann['gca_annotated'] = ann['gca']

        # Create a DataFrame with all GCA and Version pairs
        latest_versions = asm[['gca', 'version', 'gca_latest']].rename(
            columns={'version': 'latest_version'}
        )
        latest_versions['gca'] = latest_versions['gca'].str.replace(r'\.\d+$', '', regex=True)

        # Rename GCA column in filtered_df for clarity
        ann['gca'] = ann['gca'].str.replace(r'\.\d+$', '', regex=True)

        # Merge latest versions with filtered_df based on GCA
        anno_wide = ann.merge(latest_versions, on="gca", how='left')

        # Create new columns for annotated and assembly versions
        anno_wide['annotated_version'] = anno_wide['version']
        anno_wide['assembly_version'] = anno_wide['latest_version']


        # Check if the annotated GCA is the same as the latest assembly GCA
        def check_latest_annotated(row):
            if pd.isna(row['assembly_version']):
                return 'Yes, low quality assembly version'
            return 'Yes' if row['annotated_version'] == row['assembly_version'] else 'No'

        anno_wide['latest_annotated'] = anno_wide.apply(check_latest_annotated, axis=1)

        return anno_wide

    except Exception as e:
        logging.error(f"Error in check_if_gca_is_latest_annotated: {str(e)}")
    return None


def generate_report(end_date, start_date, group_name, taxon_id, bioproject_id, release_type):
    logging.info(f"Generating tables for end date: {end_date}, start date: {start_date}, group name: {group_name}, taxon id: {taxon_id}, bioproject id: {bioproject_id}, release type: {release_type}")
    df_meta_genebuild= query_meta_registry(end_date, start_date, group_name, taxon_id, bioproject_id, release_type)
    if df_meta_genebuild is None or df_meta_genebuild.empty:
        logging.error("No metadata records returned from query_meta_registry.")
        raise HTTPException(status_code=404, detail="No records found for the given filters.")

    logging.info(f"Adding additional info from beta prod server")
    anno_wide_pre_check = get_annotation_info_beta(df_meta_genebuild)

    logging.info(f"Checking if annotation is the latest GCA version")
    anno_wide = check_if_gca_is_latest_annotated(anno_wide_pre_check)

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

    # Create tables for charts
    number_of_annotations = (
        anno_wide[['gca', 'gb_status']]
        .groupby('gb_status')
        .size()
        .reset_index(name='number_of_annotations'))

    method_report = (
        anno_wide[['gca', 'annotation_method']]
        .groupby('annotation_method')
        .size()
        .reset_index(name='number_of_annotations'))

    num_unique_taxa = anno_wide['lowest_taxon_id'].nunique()
    top_3_taxa = (
        anno_wide.groupby(['scientific_name'])
        .size()
        .reset_index(name='count')
        .sort_values(by='count', ascending=False)
        .head(3)
    )

    project_report = (
        anno_wide[['gca', 'associated_project']]
        .groupby('associated_project')
        .size()
        .reset_index(name='number_of_annotations'))

    # Extract C: value
    busco_complete_series = (
        anno_wide['busco_protein']
        .str.extract(r'C:(\d+\.\d+)%')[0]
        .astype(float)
    )

    # Compute the average
    average_busco = busco_complete_series.mean()


    main_report = anno_wide[['associated_project', 'gca', 'genebuilder', 'gb_status', 'release_type', 'ftp', 'latest_annotated', 'busco_protein', 'date_completed_beta', 'release_date_beta']]


    # Transforming out of range float values that are not JSON compliant: nan
    logging.info(f"Transfroming Out of range float values that are not JSON compliant")
    anno_wide = anno_wide.apply(lambda col: col.fillna("") if col.dtype == "object" else col)
    number_of_annotations = number_of_annotations.apply(lambda col: col.fillna("") if col.dtype == "object" else col)
    method_report = method_report.apply(lambda col: col.fillna("") if col.dtype == "object" else col)
    top_3_taxa = top_3_taxa.apply(lambda col: col.fillna("") if col.dtype == "object" else col)
    project_report = project_report.apply(lambda col: col.fillna("") if col.dtype == "object" else col)
    main_report = main_report.apply(lambda col: col.fillna("") if col.dtype == "object" else col)

    return anno_wide, number_of_annotations, method_report, num_unique_taxa, top_3_taxa, project_report, average_busco, main_report