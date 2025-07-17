# app/services/assembly_service.py
import datetime
from fastapi import HTTPException
import pandas as pd
import json
import requests
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

def is_reference_genome(accession):
	"""
	Checks if a given accession is a reference genome by querying NCBI's Assembly database.

	Args:
		accession: Genome accession ID (e.g., GCF_000001405.39)

	Returns:
		True if it is a reference genome, False otherwise
	"""
	base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
	params = {
		"db": "assembly",
		"term": f"{accession}[Assembly Accession]",
		"retmode": "json"
	}
	logging.info(f"Checking if {accession} is a reference genome.")

	try:
		response = requests.get(base_url, params=params, timeout=10)
		if response.status_code != 200:
			logging.error(f"Error querying NCBI API: {response.status_code}")
			return False

		result = response.json()
		if not result["esearchresult"]["idlist"]:
			logging.warning(f"Accession {accession} not found in NCBI Assembly database.")
			return False

		# Fetch detailed assembly information
		assembly_id = result["esearchresult"]["idlist"][0]
		summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
		summary_params = {
			"db": "assembly",
			"id": assembly_id,
			"retmode": "json"
		}

		summary_response = requests.get(summary_url, params=summary_params, timeout=10)
		if summary_response.status_code != 200:
			logging.error(f"Error retrieving assembly summary: {summary_response.status_code}")
			raise HTTPException(
				status_code=400,
				detail=f"Error retrieving assembly summary: {summary_response.status_code}"
			)

		summary_data = summary_response.json()

		# Check if the assembly is labeled as a reference genome
		try:
			assembly_info = summary_data["result"][assembly_id]
			return assembly_info.get("refseq_category", "") == "reference genome"
		except Exception as e:
			logging.error(f"Unexpected error in is_reference: {e}", exc_info=True)
			raise HTTPException(
				status_code=500,
				detail=f"Internal server error occurred while processing assemblies: {str(e)}"
			)


	except HTTPException:
		# Re-raise HTTPExceptions as they are already properly formatted
		raise
	except Exception as e:
		logging.error(f"Unexpected error when checking reference genome: {e}", exc_info=True)
		raise HTTPException(
			status_code=500,
			detail=f"Internal server error occurred while processing assemblies: {str(e)}"
		)


def get_filtered_assemblies(bioproject_id, metric_thresholds, asm_level, asm_type, release_date, taxon_id,
                            current, pipeline, transc, transc_ena, non_annotated):
	"""
	Fetch all assemblies and their metrics, filter results based on given thresholds,
	and format the results with metrics as separate columns.

	Args:
		bioproject_id: List of BioProject IDs
		metric_thresholds: Dictionary of metric names and their threshold values
		all_metrics: List of all metrics to include in the results
		asm_level: List of assembly levels to filter by
		asm_type: List of assembly types to filter by
		release_date: Filter assemblies released after this date
		taxon_id: NCBI Taxon ID to filter by
		current: Whether to filter for current assemblies only
		pipeline: Which pipeline(s) to filter by
		transc: Whether to check transcriptomic data from registry
		transc_ena: Whether to check transcriptomic data from ena
		non_annotated: Only show non-annotated assemblies

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
				cursor.execute("SELECT bioproject_name FROM main_bioproject")
				known_names = {row["bioproject_name"] for row in cursor.fetchall()}

				bioproject_name = [bp for bp in bioproject_id if bp in known_names]
				bioproject_ids = [bp for bp in bioproject_id if bp not in known_names]

				if bioproject_name:
					conditions.append(f"mb.bioproject_name IN ({','.join(['%s'] * len(bioproject_name))})")
					params.extend(bioproject_name)
					logging.info(f"Filtering by BioProject names: {', '.join(bioproject_name)}")
				if bioproject_ids:
					conditions.append(f"b.bioproject_id IN ({','.join(['%s'] * len(bioproject_ids))})")
					params.extend(bioproject_ids)
					logging.info(f"Filtering by BioProject IDs: {', '.join(bioproject_ids)}")

			if release_date:
				if isinstance(release_date, pd.Timestamp):
					release_date = release_date.strftime('%Y-%m-%d')
				elif isinstance(release_date, (datetime.date, datetime.datetime)):
					release_date = release_date.strftime('%Y-%m-%d')
				conditions.append("a.release_date >= %s")
				params.append(release_date)
				logging.info(f"Filtering by release date: {release_date}")

			if current:
				conditions.append("a.is_current = 'current'")
				logging.info("Filtering for current assemblies")

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
                SELECT b.bioproject_id, mb.bioproject_name AS associated_project, a.asm_level, a.gca_chain, a.gca_version, a.asm_type, a.release_date, a.is_current,
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
                JOIN main_bioproject mb ON b.bioproject_id = mb.bioproject_id
                {where_clause}
                ORDER BY m.metrics_name;
            """

			cursor.execute(query, tuple(params))
			results = cursor.fetchall()
			# Check if we have any results from the main query
			if not results:
				raise HTTPException(
					status_code=404,
					detail="No assemblies found matching the specified criteria."
				)
			logging.info(f"Query executed successfully, retrieved {len(results)} results.")

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

		# Convert results to DataFrame and process
		df = pd.DataFrame(results)
		if df.empty:
			raise HTTPException(status_code=404, detail="No assemblies meet the given criteria.")


		# Process data
		df["release_date"] = pd.to_datetime(df["release_date"], errors="coerce")
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
		logging.debug(f"Wide results: {print(df_wide)}")

		# Convert specified metric columns to numeric
		for metric, _ in metric_thresholds.items():
			if metric in df_wide.columns:
				df_wide[metric] = pd.to_numeric(df_wide[metric], errors='coerce')
		logging.info(f"metrics converted to numeric")

		df_wide.reset_index(inplace=True)

		# Apply filtering based on thresholds
		if metric_thresholds:
			for metric, threshold in metric_thresholds.items():
				if metric in df_wide.columns:
					df_wide = df_wide[df_wide[metric] >= threshold]
					logging.info(f"Filtered {metric} >= {threshold}")

		# Filter by assembly level and type
		if asm_level:
			df_wide = df_wide[df_wide['asm_level'].isin(asm_level)]
			logging.info(f"Filtered assembly levels: {asm_level}")

		if asm_type:
			df_wide = df_wide[df_wide['asm_type'].isin(asm_type)]
			logging.info(f"Filtered assembly types: {asm_type}")

		# Check if any assemblies remain after filtering
		if df_wide.empty:
			return "No assemblies meet the given thresholds.", None, None, None, None

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
			logging.info(f"Meging transcriptomic data genus taxon id")
			df_wide = df_wide.merge(
				transcriptomic_df, left_on="genus_taxon_id", right_on="Taxon ID", how="left",
				suffixes=('_lowest', '_genus')
			)
			logging.info(f"After genus_taxon_id merge: {df_wide.shape}")
			logging.info(f"Columns: {df_wide.columns.tolist()}")

			# Drop redundant 'Taxon ID' columns (both for lowest and genus)
			df_wide.drop(columns=["Taxon ID_lowest", "Taxon ID_species", "Taxon ID"], inplace=True)

		# Filter by pipeline if requested
		if pipeline:
			logging.info(f"Filtering results by pipeline(s): {pipeline}")
			df_wide = df_wide[df_wide['pipeline'].isin(pipeline)]

			if df_wide.empty:
				raise HTTPException(
					status_code=400,
					detail=f"No assemblies matching pipeline filter: {pipeline}"
				)

		# Filter non annoteted assemblies only

		if non_annotated:
			logging.info(f"Filtring for non-annotated assemblies")
			df_wide = df_wide
			df_wide = df_wide[df_wide['annotation_status'] == 'not started']

			if df_wide.empty:
				raise HTTPException(
					status_code=400,
					detail=f"No assemblies matching pipeline filter: {non_annotated}"
				)


		df_wide = df_wide.drop_duplicates(subset='gca', keep='first')

		# Create df_main table
		df_main = df_wide[['bioproject_id', 'associated_project', 'gca', 'scientific_name', 'release_date',
		                     'lowest_taxon_id', 'genus_taxon_id', 'internal_clade', 'asm_type', 'asm_name', 'refseq_accession', 'is_current', 'asm_level',
		              'contig_n50', 'total_sequence_length']]
		df_main = df_main.drop_duplicates(subset=['gca'], keep='first')
		logging.info(f"Created main table")
		# Clean final results
		columns_to_drop = ['contig_l50', 'gc_count', 'number_of_component_sequences', 'scaffold_l50',
		                   'total_ungapped_length', 'number_of_organelles', 'total_number_of_chromosomes',
		                   'gaps_between_scaffolds_count']

		df_wide.drop(columns=columns_to_drop, inplace=True, errors='ignore')
		df_wide = df_wide.drop_duplicates(subset=['gca'], keep='first')


		# Create GCA list
		df_gca_list = df_wide[["gca"]]
		logging.info(f"Created gca_list")

		df_main = df_main.apply(lambda col: col.fillna("") if col.dtype == "object" else col)
		df_wide = df_wide.apply(lambda col: col.fillna("") if col.dtype == "object" else col)


		return df_wide, df_main, df_gca_list, taxonomy_dict


	except HTTPException:
		# Re-raise HTTPExceptions as they are already properly formatted
		raise
	except Exception as e:
		logging.error(f"Unexpected error in get_filtered_assemblies: {e}", exc_info=True)
		raise HTTPException(
			status_code=500,
			detail=f"Internal server error occurred while processing assemblies: {str(e)}"
		)