# app/services/assembly_service.py
import datetime
from http.client import HTTPException

import numpy as np
import pandas as pd
import json
import logging
from metadata_app.backend.app.core.database import get_db_connection
from metadata_app.backend.app.services.taxonomy_service import get_descendant_taxa, assign_clade_and_species, load_clade_data
from metadata_app.backend.app.services.transcriptomics_service import add_transc_data_to_df
from metadata_app.backend.app.services.get_transcriptomic_data_ENA_service import add_data_from_ena

def load_bioproject_mapping():
    """Hardcoded path for clade settings."""
    json_file = "data/bioproject_mapping.json"
    with open(json_file, "r") as f:
        logging.info("Loading bioproject mapping json file.")
        return json.load(f)


def get_filtered_assemblies(bioproject_id, candidate, taxon_id,
                            pipeline, transc, transc_ena, start_date, end_date, group_name):
    """
    Fetch all assemblies and their metrics.

    Args:
        bioproject_id: List of BioProject IDs
        candidate: Only show annotation candidate assemblies N50>100000, don't show contig level, non-annotated, current
        start_date: Filter assemblies released after this date
        end_date: Filter assemblies released before this date
        taxon_id: NCBI Taxon ID to filter by
        pipeline: Which pipeline(s) to filter by
        transc: Whether to check transcriptomic data from registry
        transc_ena: Whether to check transcriptomic data from ena
        group_name: Filter assemblies by group name

    Returns:
        df_main: DataFrame of filtered assembly metrics
        summary_df: Summary statistics DataFrame
        dr_wide: DataFrame with additional information
        df_gca_list: DataFrame with GCA accessions
        taxonomy_dict: Dictionary of taxonomy information
    """
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
                    return f"The following BioProject IDs were not found: {', '.join(invalid_bioprojects)}", None, None, None, None

            # Build query conditions
            conditions = []
            params = []

            if bioproject_id:
                conditions.append(f"b.bioproject_id IN ({','.join(['%s'] * len(bioproject_id))})")
                params.extend(bioproject_id)
                logging.info(f"Filtering by BioProject IDs: {', '.join(bioproject_id)}")

            if group_name:
                conditions.append("g.group_name = %s")
                params.append(group_name)
                logging.info(f"Filtering by group name: {group_name}")

            if start_date:
                logging.info(f"Retrieving annotation for start date {start_date}.")
                if isinstance(start_date, pd.Timestamp):
                    start_date = start_date.strftime('%Y-%m-%d')
                elif isinstance(start_date, (datetime.date, datetime.datetime)):
                    start_date = start_date.strftime('%Y-%m-%d')
                conditions.append("a.release_date >= %s")
                params.append(start_date)

            if end_date:
                logging.info(f"Retrieving annotation for end date {end_date}.")
                if isinstance(end_date, pd.Timestamp):
                    end_date = end_date.strftime('%Y-%m-%d')
                elif isinstance(end_date, (datetime.date, datetime.datetime)):
                    end_date = end_date.strftime('%Y-%m-%d')
                conditions.append("a.release_date <= %s")
                params.append(end_date)


            if taxon_id:
                all_descendant_taxa = set()
                for tax_id in taxon_id:
                    descendant_taxa = get_descendant_taxa(tax_id)
                    if not descendant_taxa:
                        logging.warning(f"No descendants found for taxon ID {tax_id}")
                    all_descendant_taxa.update(descendant_taxa)

                if not all_descendant_taxa:
                    return f"No descendant taxa found for any of the provided Taxon IDs: {', '.join(map(str, taxon_id))}", None, None, None, None

                conditions.append(f"s.lowest_taxon_id IN ({','.join(['%s'] * len(all_descendant_taxa))})")
                params.extend(all_descendant_taxa)
                logging.info(f"Filtering by lowest taxon IDs: {', '.join(str(id) for id in all_descendant_taxa)}")

            # Create WHERE clause
            where_clause = " WHERE " + " AND ".join(conditions) if conditions else ""

            # Main query
            query = f"""
                SELECT b.bioproject_id, a.asm_level, a.gca_chain, a.gca_version, a.asm_type, a.release_date, a.is_current,
                       m.metrics_name, m.metrics_value, s.scientific_name, s.common_name, a.asm_name,
                       a.lowest_taxon_id, g.group_name, a.refseq_accession, o.infra_type, o.infra_name,
                       IF(gb.assembly_id IS NOT NULL, 'annotated/started', 'not started') AS annotation_status
                FROM bioproject b
                JOIN assembly_metrics m ON b.assembly_id = m.assembly_id
                JOIN assembly a ON m.assembly_id = a.assembly_id
                LEFT JOIN species s ON a.lowest_taxon_id = s.lowest_taxon_id
                LEFT JOIN group_assembly g ON a.assembly_id = g.assembly_id
                LEFT JOIN organism o ON a.assembly_id = o.assembly_id
                LEFT JOIN genebuild gb ON a.assembly_id = gb.assembly_id
                {where_clause}
                AND a.is_current = 'current'
                ORDER BY m.metrics_name;
            """

            cursor.execute(query, tuple(params))
            results = cursor.fetchall()
            logging.info(f"Query executed successfully, retrieved {len(results)} results.")

            # Get taxonomy data
            lowest_taxon_ids = {row['lowest_taxon_id'] for row in results if
                                'lowest_taxon_id' in row and row['lowest_taxon_id'] is not None}
            logging.debug(f"Collected lowest taxon IDs {print(lowest_taxon_ids)}")

            if not lowest_taxon_ids:
                # No results or no taxon IDs found
                return "No assemblies meet the given criteria.", None, None, None, None

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

        # Convert results to DataFrame and process
        df = pd.DataFrame(results)
        if df.empty:
            return "No assemblies meet the given criteria.", None, None, None, None

        #Only show non annotated assemblies
        logging.info(f"Filtring for non-annotated assemblies")
        df = df
        df = df[df['annotation_status'] == 'not started']

        #Load bioproject_mapping
        bioproject_mapping = load_bioproject_mapping()

        # Process data
        df["release_date"] = pd.to_datetime(df["release_date"], errors="coerce")
        df["associated_project"] = df["bioproject_id"].map(bioproject_mapping)
        df["gca"] = df["gca_chain"].astype(str) + "." + df["gca_version"].astype(str)

        # Clean genome_coverage
        logging.info(f"Cleaning genome coverage")
        df['metrics_value'] = df.apply(
            lambda row: float(row['metrics_value'].rstrip('x')) if row[
                                                                       'metrics_name'] == 'genome_coverage' and isinstance(
                row['metrics_value'], str) else row['metrics_value'],
            axis=1
        )

        # Pivot data for metrics
        df["GCA"] = df["gca_chain"].astype(str) + "." + df["gca_version"].astype(str)
        df_wide = df.pivot(
            index=["bioproject_id", "associated_project", "group_name", "scientific_name", "lowest_taxon_id", "asm_level", "asm_type", "asm_name", "gca", "release_date", "refseq_accession",
                   "infra_type", "infra_name", "is_current", "annotation_status"],
            columns="metrics_name",
            values="metrics_value"
        )
        logging.debug(f"Data pivoted into {len(df_wide)} rows.")
        df_wide.reset_index(inplace=True)

        #Filter for candidate annotations if requested
        if candidate:
            df_wide['contig_n50'] = pd.to_numeric(df_wide['contig_n50'], errors='coerce')  # Convert to numbers
            df_wide = df_wide[
                (df_wide['asm_level'] != 'Contig') &
                (df_wide['contig_n50'] >= 100000)
                ]
            logging.info("Filtered for annotation candidates (n50>=100000 and asm_level!='Contig')")

        # Add clade, species, and genus information
        clade_data = load_clade_data()

        df_wide[['internal_clade', 'species_taxon_id', 'genus_taxon_id', 'pipeline']] = df_wide[
            'lowest_taxon_id'].apply(
            lambda x: pd.Series(assign_clade_and_species(x, clade_data, taxonomy_dict))
        )

        logging.info(f"Added clade data")
        logging.info(f"Changing genus id format")
        df_wide['genus_taxon_id'] = (
            pd.to_numeric(df_wide['genus_taxon_id'].replace('', pd.NA), errors='coerce')
            .astype('Int64')
        )
        logging.info(f"Changed genus id format")

        # Add transcriptomic data if needed
        if transc:
            df_wide = add_transc_data_to_df(df_wide, taxonomy_dict)

        if transc_ena:
            transcriptomic_df = add_data_from_ena(df_wide)

            # Merge for the lowest taxon ID
            logging.info(f"Merging transcriptomic data lowest taxon id")
            df_wide = df_wide.merge(
                transcriptomic_df, left_on="lowest_taxon_id", right_on="Taxon ID", how="left", suffixes=('', '_lowest')
            )
            logging.info(f"After lowest_taxon_id merge: {df_wide.shape}")
            logging.info(f"Columns: {df_wide.columns.tolist()}")

            # Merge for the species taxon ID
            logging.info(f"Merging transcriptomic data species taxon id")

            df_wide = df_wide.merge(
                transcriptomic_df, left_on="species_taxon_id", right_on="Taxon ID", how="left",
                suffixes=('_lowest', '_species')
            )
            logging.info(f"After species_taxon_id merge: {df_wide.shape}")
            logging.info(f"Columns: {df_wide.columns.tolist()}")

            # Merge for the genus taxon ID (separate column)
            logging.info(f"Merging transcriptomic data genus taxon id")
            df_wide = df_wide.merge(
                transcriptomic_df, left_on="genus_taxon_id", right_on="Taxon ID", how="left",
                suffixes=('_lowest', '_genus')
            )
            logging.info(f"After genus_taxon_id merge: {df_wide.shape}")
            logging.info(f"Columns: {df_wide.columns.tolist()}")

            # Drop redundant 'Taxon ID' columns (both for lowest and genus)
            df_wide.drop(columns=["Taxon ID_lowest", "Taxon ID_species", "Taxon ID"], inplace=True)

            #Add new summary column
            df_wide["transcriptomic_evidence"] = np.where(
                (df["Short-read paired-end illumina_lowest"] != 0) |
                ((df["Short-read paired-end illumina_lowest"] == 0) & (df["Short-read paired-end illumina"] >= 1)),
                "yes",
                "no"
            )

        # Filter by pipeline if requested
        if pipeline:
            logging.info(f"Filtering results by pipeline(s): {pipeline}")
            df_wide = df_wide[df_wide['pipeline'].isin(pipeline)]

            if df_wide.empty:
                return "No assemblies meet the pipeline filter criteria.", None, None, None, None


        #Add missing column if transc is not checked from ena
        if "transcriptomic_evidence" not in df_wide.columns:
            df_wide["transcriptomic_evidence"] = "not checked"

        # Create df_main table
        rep_asm_main = df_wide[['associated_project', 'gca', 'scientific_name', 'release_date',
                             'lowest_taxon_id', 'genus_taxon_id', 'transcriptomic_evidence', 'internal_clade', 'asm_type', 'asm_name', 'refseq_accession', 'asm_level',
                      'contig_n50', 'total_sequence_length']]
        rep_asm_main = rep_asm_main.drop_duplicates(subset=['gca'], keep='first')
        logging.info(f"Created main table")
        # Clean final results
        columns_to_drop = ['contig_l50', 'gc_count', 'number_of_component_sequences', 'scaffold_l50',
                           'total_ungapped_length', 'number_of_organelles', 'total_number_of_chromosomes',
                           'gaps_between_scaffolds_count']

        df_wide.drop(columns=columns_to_drop, inplace=True, errors='ignore')
        rep_asm_wide = df_wide.drop_duplicates(subset=['gca'], keep='first')




        return rep_asm_wide, rep_asm_main, taxonomy_dict
    except Exception as e:
        logging.error(f"Error in get_filtered_assemblies: {str(e)}")
        return f"Error processing request: {str(e)}", None, None, None, None

def generate_tables(bioproject_id, candidate, taxon_id,
                            pipeline, transc, transc_ena, start_date, end_date, group_name):
    logging.info(f"Generating tables for end date: {end_date}, start date: {start_date}, group name: {group_name}, taxon id: {taxon_id}, bioproject id: {bioproject_id}, pipeline type: {pipeline}")

    rep_asm_wide, rep_asm_main, taxonomy_dict = get_filtered_assemblies(bioproject_id, candidate, taxon_id,
                            pipeline, transc, transc_ena, start_date, end_date, group_name)
    if rep_asm_wide is None or rep_asm_wide.empty:
        logging.error("No records returned from registry.")
        raise HTTPException(status_code=404, detail="No records found for the given filters.")

    # Create GCA list
    df_gca_list = rep_asm_wide[["gca"]]
    logging.info(f"Created gca_list")

    # Create output for report
    project_report = (
        rep_asm_wide[['gca', 'associated_project']]
        .groupby('associated_project')
        .size()
        .reset_index(name='count'))

    transc_data = (
        rep_asm_wide[['gca', 'transcriptomic_evidence']]
        .groupby('transcriptomic_evidence')
        .size()
        .reset_index(name='count'))

    num_unique_taxa = rep_asm_wide['lowest_taxon_id'].nunique()

    top_3_taxa = (
        rep_asm_wide.groupby(['scientific_name'])
        .size()
        .reset_index(name='count')
        .sort_values(by='count', ascending=False)
        .head(3)
    )

    asm_type_group = (
        rep_asm_wide[['gca', 'asm_type']]
        .groupby('asm_type')
        .size()
        .reset_index(name='count'))

    clade_group = (
        rep_asm_wide[['gca', 'internal_clade']]
        .groupby('internal_clade')
        .size()
        .reset_index(name='count'))

    asm_length = (
        rep_asm_wide[['gca', 'total_sequence_length']].copy()
    )

    asm_length = asm_length.sort_values(by='total_sequence_length', ascending=False).reset_index(drop=True)
    asm_length['total_sequence_length'] = pd.to_numeric(asm_length['total_sequence_length'], errors='coerce')
    asm_length['total_sequence_length_Gb'] = asm_length['total_sequence_length'] / 1e9
    asm_length = (
        asm_length[['gca', 'total_sequence_length_Gb']]
    )
    logging.info(f"Type of length: {type(asm_length)}")
    logging.info(f"Total length: {asm_length}")



    transc_cols = [col for col in rep_asm_wide.columns if col.endswith('_transc_assess_date')]
    if not transc_cols:
        transc_reg_count = "not checked"
    else:
        total_rows = len(rep_asm_wide)
        rows_with_transc = rep_asm_wide[transc_cols].apply(
            lambda row: any(pd.notna(val) and str(val).strip() != '' for val in row), axis=1
        ).sum()
        transc_reg_count = f"{(rows_with_transc / total_rows) * 100:.1f}%"

    # Transforming out of range float values that are not JSON compliant: nan
    logging.info(f"Transfroming Out of range float values that are not JSON compliant")
    rep_asm_wide = rep_asm_wide.apply(lambda col: col.fillna("") if col.dtype == "object" else col)
    project_report = project_report.apply(lambda col: col.fillna("") if col.dtype == "object" else col)
    rep_asm_main = rep_asm_main.apply(lambda col: col.fillna("") if col.dtype == "object" else col)
    top_3_taxa = top_3_taxa.apply(lambda col: col.fillna("") if col.dtype == "object" else col)
    asm_type_group = asm_type_group.apply(lambda col: col.fillna("") if col.dtype == "object" else col)
    clade_group = clade_group.apply(lambda col: col.fillna("") if col.dtype == "object" else col)
    asm_length = asm_length.apply(lambda col: col.fillna("") if col.dtype == "object" else col)
    transc_data = transc_data.apply(lambda col: col.fillna("") if col.dtype == "object" else col)

    logging.info(f"Run successfully")

    return project_report, num_unique_taxa, transc_reg_count, top_3_taxa,asm_type_group, clade_group, asm_length, transc_data, df_gca_list, rep_asm_wide, rep_asm_main


