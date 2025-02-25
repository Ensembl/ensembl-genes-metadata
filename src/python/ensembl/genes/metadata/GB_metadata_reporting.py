import pandas as pd
import json
import pymysql
import argparse
import os
import requests

def get_descendant_taxa(taxon_id):
    """
    Retrieves all descendant taxon IDs under the given taxon ID using NCBI E-utilities.

    :param taxon_id: The parent taxon ID (e.g., 40674 for Mammalia)
    :return: A set of descendant taxon IDs
    """
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "taxonomy",
        "term": f"txid{taxon_id}[Subtree]",
        "retmode": "json",
        "retmax": 100000  # Retrieve all possible results
    }

    response = requests.get(base_url, params=params)
    if response.status_code != 200:
        print("Error retrieving taxonomic data from NCBI.")
        return set()

    try:
        result = response.json()
        taxon_ids = set(result["esearchresult"]["idlist"])
        return taxon_ids
    except KeyError:
        print("Unexpected response format from NCBI.")
        return set()


# Cache dictionary to store retrieved parent taxa
taxonomy_cache_file = "data/taxonomy_cache.json"
parent_cache = {}

def load_taxonomy_cache():
    """Load cached taxonomy data from a JSON file."""
    global parent_cache
    if os.path.exists(taxonomy_cache_file):
        with open(taxonomy_cache_file, "r") as f:
            try:
                parent_cache = json.load(f)
            except json.JSONDecodeError:
                parent_cache = {}

def save_taxonomy_cache():
    """Save taxonomy cache to a JSON file for reuse."""
    with open(taxonomy_cache_file, "w") as f:
        json.dump(parent_cache, f, indent=4)

# Load cache at script startup
load_taxonomy_cache()

# Load database credentials from external JSON file
def load_db_config():
    config_path = "conf/db_config.json"
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Database config file '{config_path}' not found.")

    with open(config_path, "r") as f:
        return json.load(f)


# Get DB credentials
db_config = load_db_config()

# Connect to DB
def connect_db():
    """Establish connection to the MySQL database."""
    return pymysql.connect(**db_config, cursorclass=pymysql.cursors.DictCursor)

# BioProject mapping
bioproject_mapping = {
    "PRJEB40665": "DToL",
    "PRJEB61747": "ERGA/BGE",
    "PRJEB43510": "ERGA",
    "PRJNA533106": "EBP",
    "PRJEB47820": "ERGA_pilot",
    "PRJEB43743": "ASG"}

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

    response = requests.get(base_url, params=params)
    if response.status_code != 200:
        print("Error querying NCBI API.")
        return False

    result = response.json()
    if not result["esearchresult"]["idlist"]:
        print("Accession not found in NCBI Assembly database.")
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
        print("Error retrieving assembly summary.")
        return False

    summary_data = summary_response.json()

    # Check if the assembly is labeled as a reference genome
    try:
        assembly_info = summary_data["result"][assembly_id]
        return assembly_info.get("refseq_category", "") == "reference genome"
    except KeyError:
        print("Unexpected response format from NCBI.")
        return False


def load_clade_data():
    """Hardcoded path for clade settings."""
    json_file = "data/clade_settings.json"
    with open(json_file, "r") as f:
        return json.load(f)


def taxonomy_api(taxon):
    """Retrieve taxonomy information from cache or NCBI API."""
    if str(taxon) in parent_cache:
        return parent_cache[str(taxon)]  # Return cached data

    # If not cached, fetch data from the API
    uri = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/{taxon}/dataset_report"

    try:
        response = requests.get(uri)
        response.raise_for_status()
        taxon_data = response.json()

        # Cache the result
        parent_cache[str(taxon)] = taxon_data
        save_taxonomy_cache()  # Save cache to file

        return taxon_data
    except requests.exceptions.RequestException as e:
        print(f"Error retrieving taxonomy data for {taxon}: {e}")
        return None


def assign_clade(lowest_taxon_id, clade_data):
    """Assign internal clade based on taxonomy. Uses cache for efficiency."""
    taxon_data = taxonomy_api(lowest_taxon_id)

    if not taxon_data:
        return "Unassigned"

    parents = taxon_data.get('reports', [{}])[0].get('taxonomy', {}).get('parents', [])

    # Check if the lowest taxon ID itself is in the clade data
    for clade_name, details in clade_data.items():
        if details.get("taxon_id") == lowest_taxon_id:
            return clade_name

    # Loop through parent taxon IDs to find a match
    for parent_taxon in parents:
        if str(parent_taxon) in parent_cache:
            # If parent is already cached, avoid extra API calls
            for clade_name, details in clade_data.items():
                if details.get("taxon_id") == parent_taxon:
                    return clade_name
        else:
            # Fetch taxonomy for parent (if not cached)
            taxonomy_api(parent_taxon)

    return "Unassigned"


def get_filtered_assemblies(bioproject_id, metric_thresholds, all_metrics, asm_level, asm_type, release_date,taxon_id):
    """Fetch all assemblies and their metrics filter results based on given thresholds,
    and format the results with metrics as separate columns."""

    conn = connect_db()
    cursor = conn.cursor()

    # Check if bioproject_id is provided before validating
    if bioproject_id:
        cursor.execute("SELECT DISTINCT bioproject_id FROM bioproject;")
        valid_bioprojects = {row['bioproject_id'] for row in cursor.fetchall()}

        # Find invalid BioProject IDs
        invalid_bioprojects = set(bioproject_id) - valid_bioprojects
        if invalid_bioprojects:
            conn.close()
            error_message = f"The following BioProject IDs were not found in the database: {', '.join(invalid_bioprojects)}"
            return error_message, None, None, None  # Return error message directly

    # Build dynamic SQL filtering
    conditions = []
    params = []

    if bioproject_id:
        conditions.append(f"b.bioproject_id IN ({','.join(['%s'] * len(bioproject_id))})")
        params.extend(bioproject_id)

    if release_date:
        conditions.append("a.release_date >= %s")
        params.append(release_date)

    if taxon_id:
        descendant_taxa = get_descendant_taxa(taxon_id)
        if not descendant_taxa:
            return f"No descendant taxa found for Taxon ID {taxon_id}.", None, None, None

        conditions.append(f"s.lowest_taxon_id IN ({','.join(['%s'] * len(descendant_taxa))})")
        params.extend(descendant_taxa)

    # If there are conditions, join them with AND; otherwise, select all
    where_clause = " WHERE " + " AND ".join(conditions) if conditions else ""

    query = f"""
        SELECT b.bioproject_id, a.asm_level, a.gca_chain, a.gca_version, a.asm_type, a.release_date, 
               m.metrics_name, m.metrics_value, s.scientific_name, s.common_name, 
               s.lowest_taxon_id, g.group_name, a.refseq_accession
        FROM bioproject b
        JOIN assembly_metrics m ON b.assembly_id = m.assembly_id
        JOIN assembly a ON m.assembly_id = a.assembly_id
        LEFT JOIN species s ON a.lowest_taxon_id = s.lowest_taxon_id
        LEFT JOIN group_assembly g ON a.assembly_id = g.assembly_id
        {where_clause}
        ORDER BY m.metrics_name;
    """
    cursor.execute(query, tuple(params))
    results = cursor.fetchall()

    # Closing connection after fetching data
    conn.close()

    # Convert results to a Pandas DataFrame
    df = pd.DataFrame(results)
    df["release_date"] = pd.to_datetime(df["release_date"], errors="coerce")


    # Add project column and GCA
    df["Associated project"] = df["bioproject_id"].map(bioproject_mapping)
    df["GCA"] = df["gca_chain"].astype(str) + "." + df["gca_version"].astype(str)

    # Clean genome_coverage by removing 'x' and converting to float
    df['metrics_value'] = df.apply(
        lambda row: float(row['metrics_value'].rstrip('x')) if row['metrics_name'] == 'genome_coverage' else row['metrics_value'], axis=1)

    # Pivot the data so each metric_name becomes a separate column and combine gca_chain and gca_version, correct date format
    df["GCA"] = df["gca_chain"].astype(str) + "." + df["gca_version"].astype(str)

    df_wide = df.pivot(index=["bioproject_id", "asm_level", "asm_type", "GCA", "release_date", "refseq_accession"], columns="metrics_name", values="metrics_value")
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


    # Apply assembly level an assembly type filtering (supporting multiple values)
    if asm_level:
        df_wide = df_wide[df_wide['asm_level'].isin(asm_level)]

    if asm_type:
        df_wide = df_wide[df_wide['asm_type'].isin(asm_type)]

    # Check if any assemblies meet the given thresholds
    if df_wide.empty:
        return "No assemblies meet the given thresholds.", None, None, None

    # Clean info results table
    df_info_result = df[['bioproject_id', 'release_date', 'scientific_name', 'common_name', 'group_name', 'Associated project', 'GCA', 'lowest_taxon_id']]
    df_info_result = df_info_result.drop_duplicates(subset=['GCA'], keep='first')

    # Drop specific columns and clean multiple GCA's
    columns_to_drop = ['contig_l50', 'release_date', 'gc_count', 'number_of_component_sequences', 'scaffold_l50', 'total_ungapped_length', 'number_of_organelles','total_number_of_chromosomes']  # Adjust this list as needed
    df_wide.drop(columns=columns_to_drop, inplace=True, errors='ignore')
    df_wide = df_wide.drop_duplicates(subset=['GCA'], keep='first')


    # Calculate summary statistics (Avg, Min, Max) for selected metrics
    summary_metrics = ['contig_n50', 'scaffold_n50', 'total_sequence_length', 'gc_percent', 'genome_coverage']
    summary_df = df_wide[summary_metrics].agg(['mean', 'min', 'max'])

    # Create the new GCA table
    df_gca_list = df_wide[["GCA"]]

    # Rename columns
    df_wide.rename(
        columns={'bioproject_id': 'BioProject ID', 'asm_level': 'Assembly level', 'number_of_contigs': 'Number of contigs', 'number_of_scaffolds': 'Number of scaffolds', 'scaffold_n50': 'Scaffold N50', 'total_sequence_length':'Sequence length', 'GCA': "GCA", 'contig_n50': 'Contig N50', 'gc_percent': 'GC%', 'genome_coverage': 'Genome coverage X', 'asm_type': "Assembly type", 'refseq_accession': "RefSeq Accession"}, inplace=True)
    summary_df.rename(columns={'genome_coverage': 'Genome coverage X', 'contig_n50': 'Conting N50', 'scaffold_n50': 'Scaffold N50',  'total_sequence_length': "Sequence length", 'gc_percent': 'GC%'}, inplace=True)
    df_info_result.rename(columns={'bioproject_id': 'BioProject ID', 'release_date': 'Release date', 'scientific_name': 'Scientific name',  'common_name': "Common name", 'group_name': 'Group name', 'lowest_taxon_id': "Lowest taxon ID"}, inplace=True)

    return df_wide, summary_df, df_info_result, df_gca_list


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
    parser.add_argument('--asm_level', type=str, nargs='+', help="Assembly level options: 'Contig', 'Scaffold', 'Chromosome', 'Complete genome'.")
    parser.add_argument('--asm_type', type=str, nargs='+', help="Assembly type: 'haploid', 'alternate-pseudohaplotype', 'unresolved-diploid', 'haploid-with-alt-loci', 'diploid'.")
    parser.add_argument('--output_dir', type=str, default='./', help="Directory to save the CSV file. Files will only be saved if assemblies meet the given thresholds.")
    parser.add_argument('--taxon_id', type=int, help="NCBI Taxon ID to filter by (e.g., 40674 for Mammalia)")

    args = parser.parse_args()

    metric_thresholds = {k: v[0] if isinstance(v, list) else v for k, v in vars(args).items() if v is not None and k not in ["bioproject_id", "output_dir", "asm_level", "asm_type", "release_date"]}

    all_metrics = ["gc_percent", "total_sequence_length", "contig_n50", "number_of_contigs", "number_of_scaffolds", "scaffold_n50", "genome_coverage"]

    df, summary_df, info_result, df_gca_list = get_filtered_assemblies(args.bioproject_id, metric_thresholds, all_metrics, args.asm_level, args.asm_type, args.release_date, args.taxon_id)


    # Check if 'df' is a string (error message)
    if isinstance(df, str):
        print(df)  # Print the error message
        return

    if df.empty:
        print("No assemblies meet the given thresholds.")
        return

    # Load clade data
    clade_data = load_clade_data()

    # Add internal clade and reference genome columns to the info_result DataFrame
    #info_result["Internal clade"] = info_result["Lowest taxon ID"].apply(lambda x: assign_clade(x, clade_data))
    info_result["Reference genome"] = info_result["GCA"].apply(is_reference_genome)

    print("Assemblies that meet the given thresholds:")
    print(df)
    print("\nSummary statistics (Avg/Min/Max) for selected metrics:")
    print(summary_df)
    print("Assembly info:")
    print(info_result)
    print("GCA list:")
    print(df_gca_list)

    # Ensure the output directory exists
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Create a base name for the files based on BioProject ID and release date
    bioproject_name = '_'.join(args.bioproject_id) if args.bioproject_id else 'all_projects'
    release_date_name = args.release_date.replace('-', '') if args.release_date else 'no_release_date'

    # Save the result as a CSV file
    output_filename = os.path.join(args.output_dir, f"{bioproject_name}_{release_date_name}_filtered_assemblies.csv")
    df.to_csv(output_filename, index=False)

    # Save the summary statistics as a CSV file
    summary_filename = os.path.join(args.output_dir, f"{bioproject_name}_{release_date_name}_filtered_assemblies_summary_statistics.csv")
    summary_df.to_csv(summary_filename, index=False)

    # Save the assembly info as a CSV file
    info_result_filename = os.path.join(args.output_dir, f"{bioproject_name}_{release_date_name}_filtered_assemblies_info_result.csv")
    info_result.to_csv(info_result_filename, index=False)

    # Save the GCA list in a file
    gca_list_filename = os.path.join(args.output_dir, f"{bioproject_name}_{release_date_name}_filtered_assemblies_gca_list.csv")
    df_gca_list.to_csv(gca_list_filename, index=False, header=False)

    print(f"\nThe results have been saved to {output_filename}")
    print(f"\nSummary statistics have been saved to {summary_filename}")
    print(f"\nAssembly info has been saved to {info_result_filename}")
    print(f"\nGCA list has been saved to {gca_list_filename}")


if __name__ == "__main__":
    main()
