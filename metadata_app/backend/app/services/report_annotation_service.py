
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
        logging.error(f"Error in check_if_gca_is_latest_annotated: {str(e)}")
        return None


def generate_report(end_date, start_date, group_name, taxon_id, bioproject_id, release_type):
    logging.info(f"Generating tables for end date: {end_date}, start date: {start_date}, group name: {group_name}, taxon id: {taxon_id}, bioproject id: {bioproject_id}, release type: {release_type}")
    df_meta_genebuild= query_meta_registry(end_date, start_date, group_name, taxon_id, bioproject_id, release_type)
    if df_meta_genebuild is None or df_meta_genebuild.empty:
        logging.error("No metadata records returned from query_meta_registry.")
        raise HTTPException(status_code=404, detail="No records found for the given filters.")

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
        df_meta_genebuild_beta_updated = pd.DataFrame(
            columns=df_meta_genebuild.columns.tolist() + ["latest_annotated"])
    logging.info(f"Updated beta: {df_meta_genebuild_beta_updated.shape}")

    # Concatenate the updated beta subset with the rest
    anno_wide_pre_check = pd.concat([df_meta_genebuild_beta_updated, df_meta_genebuild_other], ignore_index=True)
    logging.info(f"Final combined: {anno_wide_pre_check.shape}")
    logging.info(f"After beta info check: {anno_wide_pre_check.shape}")

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

    anno_wide = anno_wide.drop_duplicates(subset='gca', keep='first')

    # Create tables for charts
    # Create tables for charts
    number_of_annotations_raw = (
        anno_wide[['gca', 'gb_status']]
        .groupby('gb_status')
        .size()
        .reset_index(name='count')
    )

    # Transform into desired format
    number_of_annotations = [
        {
            "gb_status": row['gb_status'],
            "count": row['count']
        }
        for _, row in number_of_annotations_raw.iterrows()
    ]

    method_report = (
        anno_wide[['gca', 'annotation_method']]
        .groupby('annotation_method')
        .size()
        .reset_index(name='count'))

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
        .reset_index(name='count'))

    if 'busco_protein' in anno_wide.columns:
        # Extract C: value
        busco_complete_series = (
            anno_wide['busco_protein']
            .str.extract(r'C:(\d+\.\d+)%')[0]
            .astype(float)
        )

        # Compute the average
        average_busco = busco_complete_series.mean()
    else:
        anno_wide['busco_protein'] = "Not available"
        average_busco = "Not available"


    main_report = anno_wide[['associated_project', 'gca', 'genebuilder', 'gb_status', 'release_type', 'ftp', 'latest_annotated', 'busco_protein', 'date_completed_beta', 'release_date_beta']]


    # Transforming out of range float values that are not JSON compliant: nan
    logging.info(f"Transfroming Out of range float values that are not JSON compliant")
    anno_wide = anno_wide.apply(lambda col: col.fillna("") if col.dtype == "object" else col)
    method_report = method_report.apply(lambda col: col.fillna("") if col.dtype == "object" else col)
    top_3_taxa = top_3_taxa.apply(lambda col: col.fillna("") if col.dtype == "object" else col)
    project_report = project_report.apply(lambda col: col.fillna("") if col.dtype == "object" else col)
    main_report = main_report.apply(lambda col: col.fillna("") if col.dtype == "object" else col)

    logging.info(f"Run successfully")

    return anno_wide, number_of_annotations, method_report, num_unique_taxa, top_3_taxa, project_report, average_busco, main_report