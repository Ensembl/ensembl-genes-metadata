# app/services/assembly_service.py
import pandas as pd
import json
import requests
import logging
from metadata_app.backend.app.core.database import get_db_connection
from metadata_app.backend.app.services.taxonomy_service import get_descendant_taxa, assign_clade_and_species, load_clade_data
from metadata_app.backend.app.services.transcriptomics_service import add_transc_data_to_df


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
			return False

		summary_data = summary_response.json()

		# Check if the assembly is labeled as a reference genome
		try:
			assembly_info = summary_data["result"][assembly_id]
			return assembly_info.get("refseq_category", "") == "reference genome"
		except KeyError:
			logging.warning(f"Unexpected response format from NCBI.")
			return False

	except requests.exceptions.RequestException as e:
		logging.error(f"Request error when checking reference genome: {e}")
		return False
	except Exception as e:
		logging.error(f"Error checking reference genome: {e}")
		return False


def get_filtered_assemblies(bioproject_id, metric_thresholds, all_metrics, asm_level, asm_type, release_date, taxon_id,
                            current, pipeline):
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

	Returns:
		df: DataFrame of filtered assembly metrics
		summary_df: Summary statistics DataFrame
		info_result: DataFrame with additional information
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

			if release_date:
				conditions.append("a.release_date >= %s")
				params.append(release_date)
				logging.info(f"Filtering by release date: {release_date}")

			if current:
				conditions.append("a.is_current = 'current'")
				logging.info("Filtering for current assemblies")

			if taxon_id:
				descendant_taxa = get_descendant_taxa(taxon_id)
				if not descendant_taxa:
					logging.error(f"No descendant taxon IDs found for {taxon_id}.")
					return f"No descendant taxa found for Taxon ID {taxon_id}.", None, None, None, None

				conditions.append(f"s.lowest_taxon_id IN ({','.join(['%s'] * len(descendant_taxa))})")
				params.extend(descendant_taxa)
				logging.info(f"Filtering by lowest taxon ID: {', '.join(str(id) for id in descendant_taxa)}")

			# Create WHERE clause
			where_clause = " WHERE " + " AND ".join(conditions) if conditions else ""

			# Main query
			query = f"""
                SELECT b.bioproject_id, a.asm_level, a.gca_chain, a.gca_version, a.asm_type, a.release_date, a.is_current,
                       m.metrics_name, m.metrics_value, s.scientific_name, s.common_name, a.asm_name,
                       s.lowest_taxon_id, g.group_name, a.refseq_accession, o.infra_type, o.infra_name
                FROM bioproject b
                JOIN assembly_metrics m ON b.assembly_id = m.assembly_id
                JOIN assembly a ON m.assembly_id = a.assembly_id
                LEFT JOIN species s ON a.lowest_taxon_id = s.lowest_taxon_id
                LEFT JOIN group_assembly g ON a.assembly_id = g.assembly_id
                LEFT JOIN organism o ON a.assembly_id = o.assembly_id
                {where_clause}
                ORDER BY m.metrics_name;
            """

			cursor.execute(query, tuple(params))
			results = cursor.fetchall()

			# Get taxonomy data
			lowest_taxon_ids = {row['lowest_taxon_id'] for row in results if
			                    'lowest_taxon_id' in row and row['lowest_taxon_id'] is not None}

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

		#Load bioproject_mapping
		bioproject_mapping = load_bioproject_mapping()

		# Process data
		df["release_date"] = pd.to_datetime(df["release_date"], errors="coerce")
		df["associated_project"] = df["bioproject_id"].map(bioproject_mapping)
		df["gca"] = df["gca_chain"].astype(str) + "." + df["gca_version"].astype(str)

		# Clean genome_coverage
		df['metrics_value'] = df.apply(
			lambda row: float(row['metrics_value'].rstrip('x')) if row[
				                                                       'metrics_name'] == 'genome_coverage' and isinstance(
				row['metrics_value'], str) else row['metrics_value'],
			axis=1
		)

		# Pivot data for metrics
		df["GCA"] = df["gca_chain"].astype(str) + "." + df["gca_version"].astype(str)
		df_wide = df.pivot(
			index=["bioproject_id", "associated_project", "scientific_name", "lowest_taxon_id", "asm_level", "asm_type", "asm_name", "gca", "release_date", "refseq_accession",
			       "infra_type", "infra_name", "is_current"],
			columns="metrics_name",
			values="metrics_value"
		)

		# Ensure all metrics are present
		for metric in all_metrics:
			if metric not in df_wide.columns:
				df_wide[metric] = None

		# Convert metric values to numeric
		for metric in all_metrics:
			if metric in df_wide.columns:
				df_wide[metric] = pd.to_numeric(df_wide[metric], errors='coerce')

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

		# Create df_main table
		df_main = df_wide[['bioproject_id', 'associated_project', 'gca', 'scientific_name', 'release_date',
		                     'lowest_taxon_id', 'asm_type', 'asm_name', 'refseq_accession', 'is_current', 'asm_level',
		              'contig_n50', 'total_sequence_length']]
		df_main = df_main.drop_duplicates(subset=['gca'], keep='first')

		# Clean final results
		columns_to_drop = ['contig_l50', 'gc_count', 'number_of_component_sequences', 'scaffold_l50',
		                   'total_ungapped_length', 'number_of_organelles', 'total_number_of_chromosomes',
		                   'gaps_between_scaffolds_count']
		df_wide.drop(columns=columns_to_drop, inplace=True, errors='ignore')
		df_wide = df_wide.drop_duplicates(subset=['gca'], keep='first')

		# Calculate summary statistics
		summary_metrics = ['contig_n50', 'scaffold_n50', 'total_sequence_length', 'gc_percent', 'genome_coverage']
		summary_df = df_wide[summary_metrics].agg(['mean', 'min', 'max'])

		# Create GCA list
		df_gca_list = df_wide[["gca"]]

		# Filter info result to match filtered assemblies
		df_main = df_main[df_main['gca'].isin(df_wide['gca'])]

		# Add clade, species, and genus information
		clade_data = load_clade_data()
		df_main[['internal_clade', 'species_taxon_id', 'genus_taxon_id', 'pipeline']] = df_main[
			'lowest_taxon_id'].apply(
			lambda x: pd.Series(assign_clade_and_species(x, clade_data, taxonomy_dict))
		)

		df_main['genus_taxon_id'] = df_main['genus_taxon_id'].astype('Int64')  # Nullable integer type

		# Add transcriptomic data if needed
		df_main = add_transc_data_to_df(df_main, taxonomy_dict)

		# Filter by pipeline if requested
		if pipeline:
			logging.info(f"Filtering results by pipeline(s): {pipeline}")
			df_main = df_main[df_main['pipeline'].isin(pipeline)]
			df_wide = df_wide[df_wide['gca'].isin(df_main['gca'])]

			if df_main.empty:
				return "No assemblies meet the pipeline filter criteria.", None, None, None, None

		return df_wide, summary_df, df_main, df_gca_list, taxonomy_dict

	except Exception as e:
		logging.error(f"Error in get_filtered_assemblies: {str(e)}")
		return f"Error processing request: {str(e)}", None, None, None, None