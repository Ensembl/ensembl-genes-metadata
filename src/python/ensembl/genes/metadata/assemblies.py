import sys

import pandas as pd
import json
import pymysql
import argparse
import os
import requests
import asyncio
import time
import logging
from check_transcriptomic_data import check_data_from_ena

# Setup logging configuration
def setup_logging(output_dir):
    log_filename = os.path.join(output_dir, "assemblies_log.txt")
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_filename),
            logging.StreamHandler()
        ],
        force=True,
    )
    logging.info("Logging initialized")


# Load database credentials from external JSON file
def load_db_config():
    config_path = "conf/db_config.json"
    if not os.path.exists(config_path):
        logging.error("Config file not found")
        raise FileNotFoundError(f"Database config file '{config_path}' not found.")

    with open(config_path, "r") as f:
        return json.load(f)


# Get DB credentials
def connect_db(config_key):
    """Establish connection to the MySQL database."""
    db_config = load_db_config()
    return pymysql.connect(**db_config[config_key], cursorclass=pymysql.cursors.DictCursor)


# BioProject mapping
bioproject_mapping = {
    "PRJEB40665": "DToL",
    "PRJEB61747": "ERGA/BGE",
    "PRJEB43510": "ERGA",
    "PRJNA533106": "EBP",
    "PRJEB47820": "ERGA_pilot",
    "PRJEB43743": "ASG",
    "PRJNA489243": "VGP",
    "PRJNA813333": "CBP"}

def is_reference_genome(accession):
    """
    Checks if a given accession is a reference genome by querying NCBI's Assembly database.

    :param accession: Genome accession ID (e.g., GCF_000001405.39)
    :return: True if it is a reference genome, False otherwise
    """
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "assembly",
        "term": f"{accession}[Assembly Accession]",
        "retmode": "json"
    }
    logging.info(f"Checking if {accession} is a reference genome.")

    response = requests.get(base_url, params=params)
    if response.status_code != 200:
        print("Error querying NCBI API.")
        logging.error(f"Error querying NCBI API.")
        return False

    result = response.json()
    if not result["esearchresult"]["idlist"]:
        logging.warning(f"Accession not found in NCBI Assembly database.")
        return False

    # Fetch detailed assembly information
    assembly_id = result["esearchresult"]["idlist"][0]
    summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    summary_params = {
        "db": "assembly",
        "id": assembly_id,
        "retmode": "json"
    }

    summary_response = requests.get(summary_url, params=summary_params)
    if summary_response.status_code != 200:
        logging.error(f"Error retrieving assembly summary.")
        return False

    summary_data = summary_response.json()

    # Check if the assembly is labeled as a reference genome
    try:
        assembly_info = summary_data["result"][assembly_id]
        return assembly_info.get("refseq_category", "") == "reference genome"
    except KeyError:
        logging.warning(f"Unexpected response format from NCBI.")

        return False


def load_clade_data():
    """Hardcoded path for clade settings."""
    json_file = "data/clade_settings.json"
    with open(json_file, "r") as f:
        logging.info("Loading clade settings json file.")
        return json.load(f)



def assign_clade_and_species(lowest_taxon_id, clade_data, taxonomy_dict, vertebrata_taxon_id=7742, human = 9606):
    """Assign internal clade and species taxon ID based on taxonomy using the provided clade data,
       and check if the taxon ID is a descendant of the vertebrata taxon ID (7742)."""

    # Retrieve the taxonomy hierarchy from the passed dictionary (no need to query DB again)
    taxonomy_hierarchy = taxonomy_dict.get(lowest_taxon_id)

    if not taxonomy_hierarchy:
        logging.warning(f"Taxonomy hierarchy not found for taxon ID {lowest_taxon_id}")
        return "Unassigned", None, None, "anno"  # Default values

    species_taxon_id = None
    genus_taxon_id = None
    internal_clade = "Unassigned"  # Default value if no clade is found

    # Define the taxon classes in hierarchical order
    taxon_classes_order = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']

    # First pass: Set species and genus taxon IDs if available
    for taxon in taxonomy_hierarchy:
        if taxon['taxon_class'] == 'species':
            species_taxon_id = taxon['taxon_class_id']
        elif taxon['taxon_class'] == 'genus':
            genus_taxon_id = taxon['taxon_class_id']

    # Second pass: Try to assign clade
    for taxon_class in taxon_classes_order:
        matching_taxon = next((t for t in taxonomy_hierarchy if t['taxon_class'] == taxon_class), None)

        if matching_taxon:
            taxon_class_id = matching_taxon['taxon_class_id']

            # Check for matching taxon_id in clade settings
            for clade_name, details in clade_data.items():
                if details.get("taxon_id") == taxon_class_id:
                    internal_clade = clade_name
                    break

            if internal_clade != "Unassigned":
                break

    # Check if chordata is in the hierarchy
    is_chordata = any(t['taxon_class_id'] == vertebrata_taxon_id for t in taxonomy_hierarchy)
    if lowest_taxon_id == human:
        pipeline = "hprc"
    elif is_chordata:
        pipeline = "main"
    else:
        pipeline = "anno"

    # Log the assignment results for debugging
    logging.debug(
        f"Taxon {lowest_taxon_id}: clade={internal_clade}, species_id={species_taxon_id}, genus_id={genus_taxon_id}, pipeline={pipeline}")

    return internal_clade, species_taxon_id, genus_taxon_id, pipeline

def get_descendant_taxa(taxon_id):
    """
    Retrieves all descendant taxon IDs under the given taxon ID using NCBI E-utilities with pagination.

    :param taxon_id: The parent taxon ID (e.g., 40674 for Mammalia)
    :return: A set of descendant taxon IDs
    """
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "taxonomy",
        "term": f"txid{taxon_id}[Subtree]",
        "retmode": "json",
        "retmax": 100000,  # Fetch in chunks
        "retstart": 0,
        "tool": "your_tool_name",
        "email": "your_email@example.com"
    }

    taxon_ids = set()

    while True:
        response = requests.get(base_url, params=params)
        if response.status_code != 200:
            logging.error(f"Error retrieving taxonomic data from NCBI. HTTP {response.status_code}.")
            break

        try:
            result = response.json()
            batch_ids = result.get("esearchresult", {}).get("idlist", [])
            logging.info(f"Descendant taxon ID lookup successful.")
            if not batch_ids:
                break  # No more results

            taxon_ids.update(batch_ids)

            # Update retstart to fetch the next batch
            params["retstart"] += params["retmax"]

            # Respect NCBI rate limits
            time.sleep(0.5)  # Avoid overloading NCBI servers

        except Exception as e:
            print(f"Error processing response: {e}")
            logging.error(f"Error processing NCBI response: {e}")
            break

    return taxon_ids

def get_trancriptomic_assessment(taxonomy_dict):
    db_config_key_t = "transcriptomic"
    conn = connect_db(db_config_key_t)
    cursor = conn.cursor()
    logging.info(f"Retrieving trancriptomic assesment data from the registry.")

    trans_taxon_ids = set()

    for lowest_taxon_id, tax_list in taxonomy_dict.items():
        trans_taxon_ids.add(lowest_taxon_id)  # Always include lowest_taxon_id
        for taxon in tax_list:
            if taxon['taxon_class'] in ['species', 'genus']:
                trans_taxon_ids.add(taxon['taxon_class_id'])

    trans_taxon_ids = list(trans_taxon_ids)

    trans_placeholders = ','.join(['%s'] * len(trans_taxon_ids))
    where_clause = f"WHERE m.taxon_id IN ({trans_placeholders})"

    query = f"""
        SELECT m.taxon_id, m.last_check AS transc_assess_date, r.qc_status AS transc_status
        FROM meta m
        JOIN run r ON m.taxon_id = r.taxon_id
        {where_clause};
    """

    cursor.execute(query, trans_taxon_ids)
    trans_results = cursor.fetchall()

    trans_df = pd.DataFrame(trans_results, columns=['taxon_id', 'transc_assess_date', 'transc_status'])

    found_ids = trans_df['taxon_id'].nunique()
    logging.info(
        f"Transcriptomic assessment data retrieved for {found_ids} taxon IDs out of {len(trans_taxon_ids)} provided.")

    missing_count = len(trans_taxon_ids) - found_ids
    if missing_count > 0:
        logging.info(f"{missing_count} taxon IDs had no transcriptomic assessment data.")

    return trans_df


def add_transc_data_to_df(info_df, taxonomy_dict):

    trans_df = get_trancriptomic_assessment(taxonomy_dict)

    for level in ['lowest', 'species', 'genus']:
        merged = info_df[[f"{level}_taxon_id"]].merge(
            trans_df,
            how='left',
            left_on=f"{level}_taxon_id",
            right_on='taxon_id'
        )
        info_df[f"{level}_transc_assess_date"] = merged['transc_assess_date']
        info_df[f"{level}_transc_status"] = merged['transc_status']

    return info_df


def get_filtered_assemblies(bioproject_id, metric_thresholds, all_metrics, asm_level, asm_type, release_date, taxon_id, current, pipeline):
    """Fetch all assemblies and their metrics filter results based on given thresholds,
    and format the results with metrics as separate columns."""
    db_config_key = "meta"
    conn = connect_db(db_config_key)
    cursor = conn.cursor()

    # Check if bioproject_id is provided before validating
    if bioproject_id:
        cursor.execute("SELECT DISTINCT bioproject_id FROM bioproject;")
        valid_bioprojects = {row['bioproject_id'] for row in cursor.fetchall()}
        invalid_bioprojects = set(bioproject_id) - valid_bioprojects
        if invalid_bioprojects:
            conn.close()
            logging.error(f"The following BioProject IDs were not found: {', '.join(invalid_bioprojects)}")
            return f"The following BioProject IDs were not found: {', '.join(invalid_bioprojects)}", None, None, None

    # Build dynamic SQL filtering
    conditions = []
    params = []

    if bioproject_id:
        conditions.append(f"b.bioproject_id IN ({','.join(['%s'] * len(bioproject_id))})")
        params.extend(bioproject_id)
        logging.info(f"Filtering by BioProject IDs: {', '.join(bioproject_id)}")

    if release_date:
        conditions.append("a.release_date >= %s")
        params.append(release_date)
        logging.info(f"Filtering by release date: {', '.join(release_date)}")

    # Add filtering for 'current' if the argument is passed
    if current:
        conditions.append("a.is_current = 'current'")
        logging.info("Filtering by for current assemblies")

    if taxon_id:
        descendant_taxa = get_descendant_taxa(taxon_id)
        if not descendant_taxa:
            logging.error(f"No descendant taxon IDs found for {taxon_id}.")
            return f"No descendant taxa found for Taxon ID {taxon_id}.", None

        conditions.append(f"s.lowest_taxon_id IN ({','.join(['%s'] * len(descendant_taxa))})")
        params.extend(descendant_taxa)
        logging.info(f"Filtering by lowest taxon ID: {', '.join(descendant_taxa)}")


    # If there are conditions, join them with AND; otherwise, select all
    where_clause = " WHERE " + " AND ".join(conditions) if conditions else ""

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

    # Collect all unique lowest_taxon_ids from the query results
    logging.info("Fetching taxonomy table started.")
    lowest_taxon_ids = {row['lowest_taxon_id'] for row in results}

    # Fetch all taxonomy data for the collected lowest_taxon_ids at once
    taxonomy_query = """
            SELECT lowest_taxon_id, taxon_class_id, taxon_class
            FROM taxonomy
            WHERE lowest_taxon_id IN ({})
            ORDER BY FIELD(taxon_class, 'species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom');
        """.format(','.join(['%s'] * len(lowest_taxon_ids)))

    cursor.execute(taxonomy_query, tuple(lowest_taxon_ids))
    taxonomy_results = cursor.fetchall()

    # Create a dictionary to store taxonomy data based on lowest_taxon_id
    taxonomy_dict = {}
    for row in taxonomy_results:
        lowest_taxon_id = row['lowest_taxon_id']
        if lowest_taxon_id not in taxonomy_dict:
            taxonomy_dict[lowest_taxon_id] = []
        taxonomy_dict[lowest_taxon_id].append({
            'taxon_class_id': row['taxon_class_id'],
            'taxon_class': row['taxon_class']
        })
    logging.info("Fetching taxonomy table done.")

    conn.close()

    # Convert results to a Pandas DataFrame
    df = pd.DataFrame(results)
    if df.empty:
        logging.info(f"No assemblies meet the following conditions: {where_clause}.")
        return "No assemblies meet the given thresholds.", None, None, None

    # Make sure release date is in the right format
    df["release_date"] = pd.to_datetime(df["release_date"], errors="coerce")

    # Add project column and GCA
    df["associated_project"] = df["bioproject_id"].map(bioproject_mapping)
    df["gca"] = df["gca_chain"].astype(str) + "." + df["gca_version"].astype(str)

    # Clean genome_coverage by removing 'x' and converting to float
    df['metrics_value'] = df.apply(
        lambda row: float(row['metrics_value'].rstrip('x')) if row['metrics_name'] == 'genome_coverage' else row['metrics_value'], axis=1)

    # Pivot the data so each metric_name becomes a separate column and combine gca_chain and gca_version, correct date format
    df["GCA"] = df["gca_chain"].astype(str) + "." + df["gca_version"].astype(str)

    df_wide = df.pivot(index=["bioproject_id", "asm_level", "asm_type", "asm_name", "gca", "release_date", "refseq_accession", "infra_type", "infra_name", "is_current"], columns="metrics_name", values="metrics_value")
    # Ensure all requested metrics are present as columns
    for metric in all_metrics:
        if metric not in df_wide.columns:
            df_wide[metric] = None  # Fill missing metrics with NaN

    # Convert all 'metrics_value' columns to numeric (and replace non-numeric with NaN)
    for metric in all_metrics:
        if metric in df_wide.columns:
            df_wide[metric] = pd.to_numeric(df_wide[metric], errors='coerce')

    df_wide.reset_index(inplace=True)

    # Apply filtering based on user thresholds
    if metric_thresholds:
        for metric, threshold in metric_thresholds.items():
            if metric in df_wide.columns:
                df_wide = df_wide.query(f"{metric} >= {threshold}")
                logging.info(f"Filtered {metric} >= {threshold} ")


    # Apply assembly level an assembly type filtering (supporting multiple values)
    if asm_level:
        df_wide = df_wide[df_wide['asm_level'].isin(asm_level)]
        logging.info(f"Filtered {asm_level} assemblies ")

    if asm_type:
        df_wide = df_wide[df_wide['asm_type'].isin(asm_type)]
        logging.info(f"Filtered {asm_type} assemblies ")

    # Check if any assemblies meet the given thresholds
    if df_wide.empty:
        logging.info(f"No assemblies meet given thresholds.")
        return "No assemblies meet the given thresholds.", None, None, None

    # Clean info results table
    df_info_result = df[['bioproject_id', 'release_date', 'scientific_name', 'common_name', 'group_name', 'associated_project', 'gca', 'lowest_taxon_id', 'infra_type', 'infra_name', 'is_current']]
    df_info_result = df_info_result.drop_duplicates(subset=['gca'], keep='first')

    # Drop specific columns and clean multiple GCA's
    columns_to_drop = ['contig_l50', 'release_date', 'gc_count', 'number_of_component_sequences', 'scaffold_l50', 'total_ungapped_length', 'number_of_organelles','total_number_of_chromosomes', 'gaps_between_scaffolds_count', 'infra_type', 'infra_name']  # Adjust this list as needed
    df_wide.drop(columns=columns_to_drop, inplace=True, errors='ignore')
    df_wide = df_wide.drop_duplicates(subset=['gca'], keep='first')


    # Calculate summary statistics (Avg, Min, Max) for selected metrics
    summary_metrics = ['contig_n50', 'scaffold_n50', 'total_sequence_length', 'gc_percent', 'genome_coverage']
    summary_df = df_wide[summary_metrics].agg(['mean', 'min', 'max'])

    # Create the new GCA table
    df_gca_list = df_wide[["gca"]]

    # Rename columns
    df_info_result = df_info_result[df_info_result['gca'].isin(df_wide['gca'])]

    # Load clade data
    clade_data = load_clade_data()

    # Add internal clade and species taxon ID columns to the info_result DataFrame
    logging.info(f"Adding clade, pipeline, genus and species info")
    df_info_result[['internal_clade', 'species_taxon_id', 'genus_taxon_id', 'pipeline']] = df_info_result['lowest_taxon_id'].apply(lambda x: pd.Series(assign_clade_and_species(x, clade_data, taxonomy_dict)))
    df_info_result['genus_taxon_id'] = df_info_result['genus_taxon_id'].astype('Int64')  # Nullable integer type
    df_info_result = add_transc_data_to_df(df_info_result, taxonomy_dict)
    logging.debug(f"Type of genus_taxon_id column: {df_info_result['genus_taxon_id'].dtype}")
    logging.info(f"Added clade, pipeline, genus and species info")

    if pipeline:
        logging.info(f"Filtering results by pipeline(s): {pipeline}")
        df_info_result = df_info_result[df_info_result['pipeline'].isin(pipeline)]

        if df_info_result.empty:
            logging.warning("No records matched the pipeline filter. Exiting.")
            sys.exit(0)

    return df_wide, summary_df, df_info_result, df_gca_list, taxonomy_dict


def main():
    parser = argparse.ArgumentParser(description="Fetch filtered assemblies for a given BioProject.")
    parser.add_argument('--bioproject_id', type=str, nargs='+', help="One or more BioProject IDs")
    parser.add_argument('--gc_percent', type=float, help="GC percent threshold")
    parser.add_argument('--total_sequence_length', type=float, help="Total sequence length in bp")
    parser.add_argument('--contig_n50', type=float, help="Contig N50")
    parser.add_argument('--number_of_contigs', type=float, help="Number of contigs")
    parser.add_argument('--number_of_scaffolds', type=float, help="Number of scaffolds")
    parser.add_argument('--scaffold_n50', type=float, help="Scaffold N50")
    parser.add_argument('--genome_coverage', type=float, help="Genome coverage in bp")
    parser.add_argument('--release_date', type=str, help="Filter assemblies released after this date (format: YYYY-MM-DD)")
    parser.add_argument('--asm_level', type=str, nargs='+', help="Assembly level options: 'Contig' 'Scaffold' 'Chromosome' 'Complete genome'.")
    parser.add_argument('--asm_type', type=str, nargs='+', help="Assembly type: 'haploid' 'alternate-pseudohaplotype' 'unresolved-diploid' 'haploid-with-alt-loci' 'diploid'.")
    parser.add_argument('--output_dir', type=str, default='./', help="Directory to save the CSV file. Files will only be saved if assemblies meet the given thresholds.")
    parser.add_argument('--taxon_id', type=int, help="NCBI Taxon ID to filter by (e.g., 40674 for Mammalia)")
    parser.add_argument('--reference', type=int, choices=[0, 1], default=0, help="Check if GCA is a reference genome (1 for yes, 0 for no)")
    parser.add_argument('--trans', type=int, choices=[0, 1], default=0, help="Check if a taxon ID has transcriptomic data from ENA (1 for yes, 0 for no)")
    parser.add_argument('--current', type=int, choices=[0, 1], default=0, help="Check if GCA is the most current version (1 for yes, 0 for no)")
    parser.add_argument('--pipeline', type=str, choices=["anno", "main", "hprc"], nargs='+', help="Pipeline(s) to filter by: choose one or more from anno, main, hprc")

    args = parser.parse_args()

    # Ensure the output directory exists
    if not os.path.exists(args.output_dir):
        logging.info(f"Creating output directory {args.output_dir}")
        os.makedirs(args.output_dir)

    # Set up logging with the user-specified output directory
    setup_logging(args.output_dir)

    metric_thresholds = {k: v[0] if isinstance(v, list) else v for k, v in vars(args).items() if v is not None and k not in ["bioproject_id", "output_dir", "asm_level", "asm_type", "release_date", "reference", "taxon_id"]}

    all_metrics = ["gc_percent", "total_sequence_length", "contig_n50", "number_of_contigs", "number_of_scaffolds", "scaffold_n50", "genome_coverage"]

    df, summary_df, info_result, df_gca_list, taxonomy_dict = get_filtered_assemblies(args.bioproject_id, metric_thresholds, all_metrics, args.asm_level, args.asm_type, args.release_date, args.taxon_id, args.current, args.pipeline)

    logging.info(f"Filtered assemblies: {len(df)}")

    # Check for reference genome only if the user requested it
    if args.reference:
        logging.info(f"Reference genome requested")
        info_result["reference_genome"] = info_result["gca"].apply(is_reference_genome)
        logging.info(f"Reference genome check done")

    # Check for transcriptomic data only if the user requested it
    if args.trans:
        # Check transcriptomic data for each taxon_id in the dataset
        logging.info(f"Transcriptomic data check requested requested")
        # Get unique taxon IDs and filter out NaN values
        taxon_ids = [tid for tid in info_result["lowest_taxon_id"].unique() if pd.notna(tid)]
        species_taxon_ids = [tid for tid in info_result["species_taxon_id"].unique() if pd.notna(tid)]
        genus_taxon_ids = [gtid for gtid in info_result["genus_taxon_id"].unique() if pd.notna(gtid)]

        # Count NaN values and log warning if any are found
        nan_lowest_count = info_result["lowest_taxon_id"].isna().sum()
        nan_species_count = info_result["species_taxon_id"].isna().sum()
        nan_genus_count = info_result["genus_taxon_id"].isna().sum()

        if nan_lowest_count > 0:
            logging.warning(f"Found {nan_lowest_count} NA values in lowest_taxon_id column")

        if nan_lowest_count > 0:
            logging.warning(f"Found {nan_species_count} NA values in species_taxon_id column")

        if nan_genus_count > 0:
            logging.warning(f"Found {nan_genus_count} NA values in genus_taxon_id column")

        # Combine unique taxon IDs into a set and convert to integers
        all_taxon_ids = set()
        for tid in list(set(taxon_ids).union(set(genus_taxon_ids)).union(set(species_taxon_ids))):
            try:
                all_taxon_ids.add(int(tid))
            except (ValueError, TypeError) as e:
                logging.warning(f"Could not convert taxon ID {tid} to integer: {e}")

        if not all_taxon_ids:
            logging.warning("No valid taxon IDs found for transcriptomic data check")
        else:
            logging.info(f"Found {len(all_taxon_ids)} valid taxon IDs for transcriptomic data check")

        semaphore = asyncio.Semaphore(5)

        # Run transcriptomic data retrieval asynchronously
        async def fetch_transcriptomic_data():
            return await asyncio.gather(
                *[check_data_from_ena(taxon_id, tree=True, semaphore=semaphore) for taxon_id in all_taxon_ids]
            )

        transcriptomic_results = asyncio.run(fetch_transcriptomic_data())

        # Convert results into a DataFrame
        transcriptomic_df = pd.DataFrame(transcriptomic_results)

        # Merge for the lowest taxon ID
        info_result = info_result.merge(
            transcriptomic_df, left_on="lowest_taxon_id", right_on="Taxon ID", how="left", suffixes=('', '_lowest')
        )

        # Merge for the species taxon ID
        info_result = info_result.merge(
            transcriptomic_df, left_on="species_taxon_id", right_on="Taxon ID", how="left",
            suffixes=('_lowest', '_species')
        )

        # Merge for the genus taxon ID (separate column)
        info_result = info_result.merge(
            transcriptomic_df, left_on="genus_taxon_id", right_on="Taxon ID", how="left", suffixes=('_lowest', '_genus')
        )

        # Drop redundant 'Taxon ID' columns (both for lowest and genus)
        info_result.drop(columns=["Taxon ID_lowest", "Taxon ID_species", "Taxon ID"], inplace=True)


    # Check if 'df' is a string (error message)
    if isinstance(df, str):
        print(df)
        logging.error(f"Error message: {df})") # Print the error message
        return

    if df.empty:
        print("No assemblies meet the given thresholds.")
        logging.error("No assemblies met the given thresholds.")
        return

    print("Assemblies that meet the given thresholds:")
    print(df)
    print("\nSummary statistics (Avg/Min/Max) for selected metrics:")
    print(summary_df)
    print("Assembly info:")
    print(info_result)
    print("GCA list:")
    print(df_gca_list)


    # Create a base name for the files based on BioProject ID and release date
    bioproject_name = '_'.join(args.bioproject_id) if args.bioproject_id else ''
    release_date_name = args.release_date.replace('-', '') if args.release_date else ''

    # Create the output filename based on the provided values
    output_filename = os.path.join(args.output_dir,
                                   f"{bioproject_name}{'_' if bioproject_name and release_date_name else ''}{release_date_name}_filtered_assemblies.csv")
    df.to_csv(output_filename, index=False)
    logging.info(f"Wrote {output_filename}")

    # Save the summary statistics as a CSV file
    summary_filename = os.path.join(args.output_dir,
                                   f"{bioproject_name}{'_' if bioproject_name and release_date_name else ''}{release_date_name}_filtered_assemblies_summary_statistics.csv")
    summary_df.to_csv(summary_filename, index=False)
    logging.info(f"Wrote {summary_filename}")

    # Save the assembly info as a CSV file
    info_result_filename = os.path.join(args.output_dir,
                                   f"{bioproject_name}{'_' if bioproject_name and release_date_name else ''}{release_date_name}_filtered_assemblies_info_result.csv")
    info_result.to_csv(info_result_filename, index=False)
    logging.info(f"Wrote {info_result_filename}")

    # Save the GCA list in a file
    gca_list_filename = os.path.join(args.output_dir,
                                        f"{bioproject_name}{'_' if bioproject_name and release_date_name else ''}{release_date_name}_filtered_assemblies_gca_list.csv")
    df_gca_list.to_csv(gca_list_filename, index=False, header=False)
    logging.info(f"Wrote {gca_list_filename}")

    print(f"\nThe results have been saved to {output_filename}")
    print(f"\nSummary statistics have been saved to {summary_filename}")
    print(f"\nAssembly info has been saved to {info_result_filename}")
    print(f"\nGCA list has been saved to {gca_list_filename}")


if __name__ == "__main__":
    main()
