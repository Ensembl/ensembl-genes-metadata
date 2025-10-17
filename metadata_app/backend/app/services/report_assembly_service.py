# app/services/assembly_service.py
import datetime
from fastapi import HTTPException
import numpy as np
import pandas as pd
import json
import logging
from metadata_app.backend.app.core.database import get_db_connection
from metadata_app.backend.app.services.taxonomy_service import get_descendant_taxa, assign_clade_and_species, \
	load_clade_data
from metadata_app.backend.app.services.transcriptomics_service import add_transc_data_to_df
from metadata_app.backend.app.services.get_transcriptomic_data_ENA_service import add_data_from_ena


def check_dataframe_not_empty(df, description, raise_404=True):
	"""
    Check if a DataFrame is empty and raise appropriate error.

    Args:
        df: DataFrame to check
        description: Description of the DataFrame for error message
        raise_404: If True, raise HTTPException with 404, otherwise log error and return False

    Returns:
        bool: True if DataFrame is not empty

    Raises:
        HTTPException: If DataFrame is empty and raise_404 is True
    """
	if df.empty:
		error_msg = f"DataFrame is empty: {description}"
		logging.error(error_msg)
		if raise_404:
			raise HTTPException(
				status_code=404,
				detail=f"No data found: {description}"
			)
		return False

	logging.info(f"DataFrame validated - {description}: {len(df)} rows")
	return True


def get_filtered_assemblies(bioproject_id, candidate, taxon_id,
                             transc, transc_ena, start_date, end_date, group_name):
	"""
    Fetch all assemblies and their metrics.

    Args:
        bioproject_id: List of BioProject IDs
        candidate: Only show annotation candidate assemblies N50>100000, don't show contig level, non-annotated, current
        start_date: Filter assemblies released after this date
        end_date: Filter assemblies released before this date
        taxon_id: NCBI Taxon ID to filter by
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
				conditions.append("a.release_date <= %s")
				params.append(start_date)

			if end_date:
				logging.info(f"Retrieving annotation for end date {end_date}.")
				if isinstance(end_date, pd.Timestamp):
					end_date = end_date.strftime('%Y-%m-%d')
				elif isinstance(end_date, (datetime.date, datetime.datetime)):
					end_date = end_date.strftime('%Y-%m-%d')
				conditions.append("a.release_date >= %s")
				params.append(end_date)

			if taxon_id:
				all_descendant_taxa = set()
				for tax_id in taxon_id:
					descendant_taxa = get_descendant_taxa(tax_id)
					if not descendant_taxa:
						raise ValueError(f"No descendant taxa found for Taxon ID {taxon_id}.")
					all_descendant_taxa.update(descendant_taxa)

				if not all_descendant_taxa:
					raise ValueError(
						f"No descendant taxa found for any of the provided Taxon IDs: {', '.join(map(str, taxon_id))}")

				conditions.append(f"s.lowest_taxon_id IN ({','.join(['%s'] * len(all_descendant_taxa))})")
				params.extend(all_descendant_taxa)
				logging.info(f"Filtering by lowest taxon IDs: {', '.join(str(id) for id in all_descendant_taxa)}")

				if not all_descendant_taxa:
					raise HTTPException(
						status_code=404,
						detail=f"No descendant taxa found for any of the provided Taxon IDs: {', '.join(map(str, taxon_id))}"
					)

				conditions.append(f"s.lowest_taxon_id IN ({','.join(['%s'] * len(all_descendant_taxa))})")
				params.extend(all_descendant_taxa)
				logging.info(f"Filtering by lowest taxon IDs: {', '.join(str(id) for id in all_descendant_taxa)}")

			# Create WHERE clause
			where_clause = " WHERE " + " AND ".join(conditions) if conditions else ""

			# Main query
			query = f"""
                SELECT b.bioproject_id, mb.bioproject_name AS associated_project, a.asm_level, a.gca_chain, a.gca_version, a.asm_type, a.release_date, a.is_current,
                       m.metrics_name, m.metrics_value, s.scientific_name, s.common_name, a.asm_name,
                       a.lowest_taxon_id, g.group_name, a.refseq_accession, o.infra_type, o.infra_name
                FROM bioproject b
                JOIN main_bioproject mb ON mb.bioproject_id = b.bioproject_id
                JOIN assembly_metrics m ON b.assembly_id = m.assembly_id
                JOIN assembly a ON m.assembly_id = a.assembly_id
                LEFT JOIN species s ON a.lowest_taxon_id = s.lowest_taxon_id
                LEFT JOIN custom_group g
				  ON (
				       (g.group_type = 'taxon' AND a.lowest_taxon_id = g.item)
				       OR
				       (g.group_type = 'assembly' AND a.gca_chain = g.item)
				     )
                LEFT JOIN organism o ON a.assembly_id = o.assembly_id
                LEFT JOIN genebuild_status gb ON a.assembly_id = gb.assembly_id
                {where_clause}
                AND a.is_current = 'current'
                AND gb.genebuild_status_id IS NULL
                ORDER BY m.metrics_name;
            """

			cursor.execute(query, tuple(params))
			results = cursor.fetchall()
			logging.info(f"Query executed successfully, retrieved {len(results)} results.")

			# Check if we have any results from the main query
			if not results:
				raise HTTPException(
					status_code=404,
					detail="No assemblies found matching the specified criteria."
				)

			# Get taxonomy data
			lowest_taxon_ids = {row['lowest_taxon_id'] for row in results if
			                    'lowest_taxon_id' in row and row['lowest_taxon_id'] is not None}
			logging.debug(f"Collected lowest taxon IDs: {lowest_taxon_ids}")

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

		# Convert results to DataFrame and process
		df = pd.DataFrame(results)
		check_dataframe_not_empty(df, "initial query results")

		# Process data
		df["release_date"] = pd.to_datetime(df["release_date"], errors="coerce")
		df["gca"] = df["gca_chain"].astype(str) + "." + df["gca_version"].astype(str)

		# Clean genome_coverage
		logging.info("Cleaning genome coverage")
		df['metrics_value'] = df.apply(
			lambda row: float(row['metrics_value'].rstrip('x')) if row[
				                                                       'metrics_name'] == 'genome_coverage' and isinstance(
				row['metrics_value'], str) else row['metrics_value'],
			axis=1
		)

		# Pivot data for metrics
		df["GCA"] = df["gca_chain"].astype(str) + "." + df["gca_version"].astype(str)
		df_wide = df.pivot(
			index=["bioproject_id", "associated_project", "group_name", "scientific_name", "lowest_taxon_id",
			       "asm_level", "asm_type", "asm_name", "gca", "release_date", "refseq_accession",
			       "infra_type", "infra_name", "is_current"],
			columns="metrics_name",
			values="metrics_value"
		)
		logging.debug(f"Data pivoted into {len(df_wide)} rows.")
		df_wide.reset_index(inplace=True)
		check_dataframe_not_empty(df_wide, "data after pivoting")

		# Filter for candidate annotations if requested
		if candidate:
			df_wide['contig_n50'] = pd.to_numeric(df_wide['contig_n50'], errors='coerce')  # Convert to numbers
			df_wide = df_wide[
				(~df_wide['asm_level'].isin(['Contig', 'Scaffold'])) &
				(df_wide['contig_n50'] >= 100000)
				]
			logging.info("Filtered for annotation candidates (n50>=100000 and asm_level!='Contig')")
			check_dataframe_not_empty(df_wide, "candidate assemblies (N50 >= 100000, non-contig)")

		# Add clade, species, and genus information
		clade_data = load_clade_data()

		df_wide[['internal_clade', 'species_taxon_id', 'genus_taxon_id']] = df_wide[
			'lowest_taxon_id'].apply(
			lambda x: pd.Series(assign_clade_and_species(x, clade_data, taxonomy_dict))
		)

		logging.info("Added clade data")
		logging.info("Changing genus id format")
		df_wide['genus_taxon_id'] = (
			pd.to_numeric(df_wide['genus_taxon_id'].replace('', pd.NA), errors='coerce')
			.astype('Int64')
		)
		logging.info("Changed genus id format")

		# Add transcriptomic data if needed
		if transc:
			logging.info("Adding transcriptomic data from registry")
			df_wide = add_transc_data_to_df(df_wide, taxonomy_dict)
			check_dataframe_not_empty(df_wide, "data after adding transcriptomic data from registry")

		if transc_ena:
			logging.info("Adding transcriptomic data from ENA")
			transcriptomic_df = add_data_from_ena(df_wide)

			if transcriptomic_df is not None:
				check_dataframe_not_empty(transcriptomic_df, "transcriptomic data from ENA")

				# Merge for the lowest taxon ID
				logging.info("Merging transcriptomic data lowest taxon id")
				df_wide = df_wide.merge(
					transcriptomic_df, left_on="lowest_taxon_id", right_on="Taxon ID", how="left",
					suffixes=('', '_lowest')
				)
				logging.info(f"After lowest_taxon_id merge: {df_wide.shape}")
				check_dataframe_not_empty(df_wide, "data after merging ENA transcriptomic data (lowest taxon)")

				# Merge for the species taxon ID
				logging.info("Merging transcriptomic data species taxon id")
				df_wide = df_wide.merge(
					transcriptomic_df, left_on="species_taxon_id", right_on="Taxon ID", how="left",
					suffixes=('_lowest', '_species')
				)
				logging.info(f"After species_taxon_id merge: {df_wide.shape}")
				check_dataframe_not_empty(df_wide, "data after merging ENA transcriptomic data (species taxon)")

				# Merge for the genus taxon ID (separate column)
				logging.info("Merging transcriptomic data genus taxon id")
				df_wide = df_wide.merge(
					transcriptomic_df, left_on="genus_taxon_id", right_on="Taxon ID", how="left",
					suffixes=('_lowest', '_genus')
				)
				logging.info(f"After genus_taxon_id merge: {df_wide.shape}")
				check_dataframe_not_empty(df_wide, "data after merging ENA transcriptomic data (genus taxon)")

				# Drop redundant 'Taxon ID' columns (both for lowest and genus)
				df_wide.drop(columns=["Taxon ID_lowest", "Taxon ID_species", "Taxon ID"], inplace=True)

				# Add new summary column
				df_wide["transcriptomic_evidence"] = np.where(
					(df_wide["Short-read paired-end illumina_lowest"] != 0) |
					((df_wide["Short-read paired-end illumina_lowest"] == 0) & (
								df_wide["Short-read paired-end illumina"] >= 1)),
					"yes",
					"no"
				)
			else:
				logging.warning("No transcriptomic data retrieved from ENA")


		# Add missing column if transc is not checked from ena
		if "transcriptomic_evidence" not in df_wide.columns:
			df_wide["transcriptomic_evidence"] = "not checked"


		# Clean final results
		columns_to_drop = ['contig_l50', 'gc_count', 'number_of_component_sequences', 'scaffold_l50',
		                   'total_ungapped_length', 'number_of_organelles', 'total_number_of_chromosomes',
		                   'gaps_between_scaffolds_count']

		df_wide.drop(columns=columns_to_drop, inplace=True, errors='ignore')
		agg_dict = {
			'bioproject_id': lambda x: ', '.join(sorted(set(x))),
			'associated_project': lambda x: ', '.join(sorted(set(x))),
		}

		# add all other columns except gca, bioproject_id, associated_project with 'first'
		for col in df_wide.columns:
			if col not in ['gca', 'bioproject_id', 'associated_project']:
				agg_dict[col] = 'first'

		rep_asm_wide = df_wide.groupby('gca', as_index=False).agg(agg_dict)
		check_dataframe_not_empty(rep_asm_wide, "wide results table after deduplication and cleanup")

		# Create df_main table
		rep_asm_main = df_wide[['associated_project', 'gca', 'scientific_name', 'release_date',
		                        'lowest_taxon_id', 'genus_taxon_id', 'transcriptomic_evidence', 'internal_clade',
		                        'asm_type', 'asm_name', 'refseq_accession', 'asm_level',
		                        'contig_n50', 'total_sequence_length']]
		check_dataframe_not_empty(rep_asm_main, "main results table after deduplication")
		logging.info("Created main table")

		return rep_asm_wide, rep_asm_main, taxonomy_dict

	except HTTPException:
		# Re-raise HTTPExceptions as they are already properly formatted
		raise
	except Exception as e:
		logging.error(f"Unexpected error in get_filtered_assemblies: {e}", exc_info=True)
		raise HTTPException(
			status_code=500,
			detail=f"Internal server error occurred while processing assemblies: {str(e)}"
		)


def generate_tables(bioproject_id, candidate, taxon_id,
                    transc, transc_ena, start_date, end_date, group_name):
	logging.info(
		f"Generating tables for end date: {end_date}, start date: {start_date}, group name: {group_name}, taxon id: {taxon_id}, bioproject id: {bioproject_id}")

	try:
		rep_asm_wide, rep_asm_main, taxonomy_dict = get_filtered_assemblies(bioproject_id, candidate, taxon_id,
		                                                                    transc, transc_ena, start_date,
		                                                                    end_date, group_name)
	except HTTPException:
		logging.error("HTTPException raised during assembly filtering")
		raise
	except Exception as e:
		logging.error("Unexpected error occurred during assembly filtering", exc_info=True)
		raise HTTPException(status_code=500, detail="An unexpected error occurred.")

	# Validate the returned DataFrames
	check_dataframe_not_empty(rep_asm_wide, "wide assembly results from filtering")
	check_dataframe_not_empty(rep_asm_main, "main assembly results from filtering")

	# Create GCA list
	df_gca_list = rep_asm_wide[["gca"]]
	check_dataframe_not_empty(df_gca_list, "GCA list")
	logging.info("Created gca_list")

	# Create output for report
	project_report = (
		rep_asm_wide[['gca', 'associated_project']]
		.groupby('associated_project')
		.size()
		.reset_index(name='count'))
	check_dataframe_not_empty(project_report, "project report summary")

	transc_data = (
		rep_asm_wide[['gca', 'transcriptomic_evidence']]
		.groupby('transcriptomic_evidence')
		.size()
		.reset_index(name='count'))
	check_dataframe_not_empty(transc_data, "transcriptomic data summary")

	num_unique_taxa = rep_asm_wide['lowest_taxon_id'].nunique()
	if num_unique_taxa == 0:
		logging.error("No unique taxa found")
		raise HTTPException(status_code=404, detail="No unique taxa found in the results")

	top_3_taxa = (
		rep_asm_wide.groupby(['scientific_name'])
		.size()
		.reset_index(name='count')
		.sort_values(by='count', ascending=False)
		.head(3)
	)
	check_dataframe_not_empty(top_3_taxa, "top 3 taxa summary")

	asm_type_group = (
		rep_asm_wide[['gca', 'asm_type']]
		.groupby('asm_type')
		.size()
		.reset_index(name='count'))
	check_dataframe_not_empty(asm_type_group, "assembly type summary")

	clade_group = (
		rep_asm_wide[['gca', 'internal_clade']]
		.groupby('internal_clade')
		.size()
		.reset_index(name='count'))
	check_dataframe_not_empty(clade_group, "clade group summary")

	asm_length = (
		rep_asm_wide[['gca', 'total_sequence_length']].copy()
	)
	check_dataframe_not_empty(asm_length, "assembly length data")

	asm_length = asm_length.sort_values(by='total_sequence_length', ascending=False).reset_index(drop=True)
	asm_length['total_sequence_length'] = pd.to_numeric(asm_length['total_sequence_length'], errors='coerce')
	asm_length['total_sequence_length_Gb'] = asm_length['total_sequence_length'] / 1e9
	asm_length = (
		asm_length[['gca', 'total_sequence_length_Gb']]
	)
	check_dataframe_not_empty(asm_length, "processed assembly length data")
	logging.info(f"Assembly length data processed successfully: {len(asm_length)} rows")

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
	logging.info("Transforming out of range float values that are not JSON compliant")
	rep_asm_wide = rep_asm_wide.apply(lambda col: col.fillna("") if col.dtype == "object" else col)
	project_report = project_report.apply(lambda col: col.fillna("") if col.dtype == "object" else col)
	rep_asm_main = rep_asm_main.apply(lambda col: col.fillna("") if col.dtype == "object" else col)
	top_3_taxa = top_3_taxa.apply(lambda col: col.fillna("") if col.dtype == "object" else col)
	asm_type_group = asm_type_group.apply(lambda col: col.fillna("") if col.dtype == "object" else col)
	clade_group = clade_group.apply(lambda col: col.fillna("") if col.dtype == "object" else col)
	asm_length = asm_length.apply(lambda col: col.fillna("") if col.dtype == "object" else col)
	transc_data = transc_data.apply(lambda col: col.fillna("") if col.dtype == "object" else col)

	# Final validation of all return DataFrames
	check_dataframe_not_empty(project_report, "final project report")
	check_dataframe_not_empty(top_3_taxa, "final top 3 taxa")
	check_dataframe_not_empty(asm_type_group, "final assembly type group")
	check_dataframe_not_empty(clade_group, "final clade group")
	check_dataframe_not_empty(asm_length, "final assembly length")
	check_dataframe_not_empty(transc_data, "final transcriptomic data")
	check_dataframe_not_empty(df_gca_list, "final GCA list")
	check_dataframe_not_empty(rep_asm_wide, "final wide assembly results")
	check_dataframe_not_empty(rep_asm_main, "final main assembly results")

	logging.info("All tables generated successfully")

	return project_report, num_unique_taxa, transc_reg_count, top_3_taxa, asm_type_group, clade_group, asm_length, transc_data, df_gca_list, rep_asm_wide, rep_asm_main

def generate_overview():
    try:
        with get_db_connection("meta") as conn:
            cursor = conn.cursor()

            # === 1. MAIN QUERY ===
            query = """
                SELECT 
                    COALESCE(g.group_name, mb.bioproject_name) AS project_name,
                    a.assembly_id,
                    a.asm_level,
                    m.metrics_name,
                    m.metrics_value,
                    a.lowest_taxon_id,
                    s.species_taxon_id,
                    a.is_current,
                    t.taxon_class_id AS genus_taxon_id,
                    gb.gb_status AS genebuild_status
                FROM assembly a
                LEFT JOIN bioproject b ON a.assembly_id = b.assembly_id
                LEFT JOIN main_bioproject mb ON mb.bioproject_id = b.bioproject_id
                LEFT JOIN assembly_metrics m ON b.assembly_id = m.assembly_id
                JOIN taxonomy t ON a.lowest_taxon_id = t.lowest_taxon_id
                JOIN species s ON a.lowest_taxon_id = s.lowest_taxon_id
                LEFT JOIN custom_group g
                    ON (
                        (g.group_type = 'taxon' AND a.lowest_taxon_id = g.item)
                        OR
                        (g.group_type = 'assembly' AND a.gca_chain = g.item)
                    )
                LEFT JOIN genebuild_status gb ON a.assembly_id = gb.assembly_id
                WHERE (mb.bioproject_id IS NOT NULL OR g.group_name IS NOT NULL)
                  AND t.taxon_class = "genus"
                  AND m.metrics_name IN ('contig_n50');
            """

            cursor.execute(query)
            results = cursor.fetchall()

        if not results:
            raise HTTPException(status_code=404, detail="No assemblies found matching criteria.")

        df = pd.DataFrame(results)
        logging.info(f"Fetched {len(df)} rows from database")
        print(f"Projects before pivot: {df['project_name'].nunique()}")

        df["genebuild_status"] = df["genebuild_status"].fillna("not_annotated")
        df = df.drop_duplicates(
	        subset=["project_name", "assembly_id", "lowest_taxon_id"],
	        keep="first"
        )

        # === 2. Pivot metrics wide ===
        df_wide = df.pivot_table(
            index=["project_name", "assembly_id", "asm_level", "is_current",
                   "lowest_taxon_id", "species_taxon_id", "genus_taxon_id", "genebuild_status"],
            columns="metrics_name",
            values="metrics_value",
            aggfunc="first",
	        fill_value=np.nan
        ).reset_index()

        print(f"Projects after pivot: {df_wide['project_name'].nunique()}")

        # === 3. Add ENA transcriptomic data ===
        transcriptomic_df = add_data_from_ena(df_wide)
        if transcriptomic_df is not None and not transcriptomic_df.empty:
            df_wide = df_wide.merge(
                transcriptomic_df[["Taxon ID", "Short-read paired-end illumina"]],
                left_on="genus_taxon_id", right_on="Taxon ID", how="left"
            )
            df_wide["transcriptomic_evidence"] = np.where(
                df_wide["Short-read paired-end illumina"].fillna(0) > 0, "yes", "no"
            )
        else:
            df_wide["transcriptomic_evidence"] = "no"

        # === 4. Compute per-project summary ===
        df_wide["is_annotation_candidate"] = (
            df_wide["asm_level"].isin(["Complete genome", "Chromosome"]) &
            (df_wide["contig_n50"].astype(float) > 100000) &
            (df_wide["transcriptomic_evidence"] == "yes") &
            (df_wide["is_current"] == "current")
        )

        df_wide["unannotated"] = (
		        df_wide["asm_level"].isin(["Complete genome", "Chromosome"]) &
		        (df_wide["contig_n50"].astype(float) > 100000) &
		        (df_wide["transcriptomic_evidence"] == "yes") &
		        (df_wide["is_current"] == "current") &
		        (df_wide["genebuild_status"] == "not_annotated")
        )


        # Count totals per project
        summary = (
            df_wide.groupby("project_name", as_index=False)
            .agg(
                total_assemblies=("assembly_id", "nunique"),
                annotation_candidates=("is_annotation_candidate", "sum"),
	            unannotated=("unannotated", "sum"),
                in_progress=("genebuild_status", lambda x: x.isin(["in_progress", "complete", "pre_released"]).sum()),
                live=("genebuild_status", lambda x: (x == "live").sum()),
            )
        )

        logging.info(f"Generated overview summary for {len(summary)} projects")
        summary = summary.to_dict(orient="records")
        return summary

    except Exception as e:
        logging.exception("Error generating overview table")
        raise HTTPException(status_code=500, detail=str(e))


def project_per_year():
	try:
		with get_db_connection("meta") as conn:
			cursor = conn.cursor()

			query = """
                SELECT 
                    COALESCE(g.group_name, mb.bioproject_name) AS project_name,
                    a.assembly_id,
                    a.release_date,
                    gb.release_date AS release_date_anno,
                    gb.gb_status, 
                    a.is_current
                FROM assembly a
                LEFT JOIN bioproject b ON b.assembly_id = a.assembly_id
                LEFT JOIN main_bioproject mb ON mb.bioproject_id = b.bioproject_id
                LEFT JOIN custom_group g
                    ON (
                        (g.group_type = 'taxon' AND a.lowest_taxon_id = g.item)
                        OR
                        (g.group_type = 'assembly' AND a.gca_chain = g.item)
                    )
                LEFT JOIN genebuild_status gb ON a.assembly_id = gb.assembly_id
                WHERE (mb.bioproject_id IS NOT NULL OR g.group_name IS NOT NULL)
            """

			cursor.execute(query)
			results = cursor.fetchall()

		if not results:
			raise HTTPException(status_code=404, detail="No assemblies found matching criteria.")

		df = pd.DataFrame(results)
		logging.info(f"Fetched {len(df)} rows from database")
		print(f"Projects per year before: {df['project_name'].nunique()}")


		df["gb_status"] = df["gb_status"].fillna("not_annotated")
		df = df.drop_duplicates(
			subset=["project_name", "assembly_id"],
			keep="first"
		)
		print(f"Projects per year after drop duplicates: {df['project_name'].nunique()}")

		# Ensure release_date columns are datetime
		df['release_date'] = pd.to_datetime(df['release_date'], errors='coerce')
		df['release_date_anno'] = pd.to_datetime(df['release_date_anno'], errors='coerce')
		print(f"Projects per year after datetime: {df['project_name'].nunique()}")

		# Extract year
		df_assembly=df
		df_assembly['release_year'] = df_assembly['release_date'].dt.year
		df_assembly = df_assembly[df_assembly["release_year"] >= 2019]

		df_live=df
		# Extract annotation release year, keep NaT as a placeholder
		df_live['release_year_anno'] = df_live['release_date_anno'].dt.year

		print(f"Projects per year after extract year: {df_live['project_name'].nunique()}")
		print(df_live)
		# --- 1. Total assemblies per project per year ---
		df_assembly = df_assembly[df_assembly["is_current"] == "current"]
		df_assembly = (
			df_assembly.groupby(['project_name', 'release_year'], as_index=False)
			.agg(total_assemblies=('assembly_id', 'nunique'))
		)

		# --- 2. Live annotations per project per year ---
		df_live = df_live[df_live['gb_status'] == 'live'].copy()
		print(f"Projects per year after live filter by annotations: {df_live['project_name'].nunique()}")

		df_annotations = (
			df_live.groupby(['project_name', 'release_year_anno'], as_index=False)
			.agg(live_annotations=('assembly_id', 'nunique'))
			.rename(columns={'release_year_anno': 'release_year'})
		)
		print(f"Projects per year after group by annotations: {df_annotations['project_name'].nunique()}")


		# Pivot assemblies for grouped bar chart
		df_assembly = df_assembly.pivot_table(
			index='release_year',
			columns='project_name',
			values='total_assemblies',
			fill_value=0
		).reset_index()
		df_assembly.columns.name = None  # remove columns name

		# Pivot annotations for grouped bar chart
		df_annotations = df_annotations.pivot_table(
			index='release_year',
			columns='project_name',
			values='live_annotations',
			fill_value=0
		).reset_index()
		df_annotations.columns.name = None

		df_assembly['release_year'] = df_assembly['release_year'].astype(int)
		df_annotations['release_year'] = df_annotations['release_year'].astype(int)

		df_assembly = df_assembly.to_dict(orient="records")
		df_annotations = df_annotations.to_dict(orient="records")
		print("Assemblies bar")
		print(df_assembly)
		print("Annotations bar")
		print(df_annotations)


		return df_assembly, df_annotations

	except Exception as e:
		logging.exception("Error generating project per year tables")
		raise HTTPException(status_code=500, detail=str(e))



