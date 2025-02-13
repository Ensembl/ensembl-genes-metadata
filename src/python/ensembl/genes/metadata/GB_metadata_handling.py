import pandas as pd
import json
import pymysql
import argparse
import os
import requests
from datetime import datetime


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
    "PRJEB40665": "dtol",
    "PRJEB61747": "ERGA/BGE",
    "PRJEB43510": "ERGA",
    "PRJNA533106": "EBP",}

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
    uri = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/{taxon}/dataset_report"
    response = requests.get(uri)
    response.raise_for_status()
    return response.json()

# Assign internal clade to GCA
def assign_clade(lowest_taxon_id, clade_data):
    """Assign internal clade to GCAs based on JSON file. Function looks for the taxonomy id associated to the GCAs
    then loops though parents until it finds a match in the JSON file."""

    taxon_data = taxonomy_api(lowest_taxon_id)
    parents = taxon_data['reports'][0]['taxonomy']['parents']

    # Check lowest taxon ID in the JSON data
    for clade_name, details in clade_data.items():
        if details.get("taxon_id") == lowest_taxon_id:
            return clade_name

    # Loop through parent taxon IDs to find a match
    for parent_taxon in parents:
        for clade_name, details in clade_data.items():
            if details.get("taxon_id") == parent_taxon:
                return clade_name
    return "Unassigned"

def get_filtered_assemblies(bioproject_id, metric_thresholds, all_metrics, asm_level, asm_type, release_date):
    """Fetch all assemblies and their metrics for a given BioProject ID, filter results based on given thresholds,
    and format the results with metrics as separate columns."""

    conn = connect_db()
    cursor = conn.cursor()

    # Step 1: Retrieve valid BioProject IDs from the database
    cursor.execute("SELECT DISTINCT bioproject_id FROM bioproject;")
    valid_bioprojects = {row['bioproject_id'] for row in cursor.fetchall()}

    # Step 2: Check if all user-provided BioProject IDs exist
    invalid_bioprojects = set(bioproject_id) - valid_bioprojects
    if invalid_bioprojects:
        print(f"The following BioProject IDs were not found in the database: {', '.join(invalid_bioprojects)}")
        conn.close()
        exit()

    # Step 3: Proceed with fetching and filtering the data
    # Query to fetch all metrics for the given BioProject ID
    query = f"""
            SELECT b.bioproject_id, a.asm_level, a.gca_chain, a.gca_version, a.asm_type, a.release_date, m.metrics_name, m.metrics_value
            FROM bioproject b
            JOIN assembly_metrics m ON b.assembly_id = m.assembly_id
            JOIN assembly a ON m.assembly_id = a.assembly_id
            WHERE b.bioproject_id IN ({",".join(["%s"] * len(bioproject_id))})
            ORDER BY m.metrics_name;
        """

    # Query to fetch refseq for the given BioProject ID
    refseq_count_query = """
        SELECT COUNT(*) AS refseq_count
        FROM assembly a
        JOIN bioproject b ON a.assembly_id = b.assembly_id
        WHERE b.bioproject_id IN ({}) AND a.refseq_accession IS NOT NULL AND a.refseq_accession != '';
    """

    # Query to fetch all extra information for the given BioProject ID
    info_query = """
        SELECT b.bioproject_id, s.scientific_name, s.common_name, s.lowest_taxon_id, g.group_name
        FROM species s
        LEFT JOIN assembly a ON s.lowest_taxon_id = a.lowest_taxon_id
        LEFT JOIN group_assembly g ON a.assembly_id = g.assembly_id
        LEFT JOIN bioproject b ON a.assembly_id = b.assembly_id
        WHERE b.bioproject_id IN ({})
        ORDER BY s.scientific_name;
    """
    query = query.format(",".join(["%s"] * len(bioproject_id)))
    cursor.execute(query, tuple(bioproject_id))
    results = cursor.fetchall()

    refseq_count_query = refseq_count_query.format(",".join(["%s"] * len(bioproject_id)))
    cursor.execute(refseq_count_query, tuple(bioproject_id))
    refseq_count_result = cursor.fetchone()
    refseq_count = refseq_count_result['refseq_count']

    info_query = info_query.format(",".join(["%s"] * len(bioproject_id)))
    cursor.execute(info_query, tuple(bioproject_id))
    info_result = cursor.fetchall()

    conn.close()

    # Convert results to a Pandas DataFrame
    df = pd.DataFrame(results)
    df["release_date"] = pd.to_datetime(df["release_date"], errors="coerce")
    df_info_result = pd.DataFrame(info_result)

    # Add project column and GCA ID to info table based on BioProject_ID
    df_info_result["Project"] = df_info_result["bioproject_id"].map(bioproject_mapping)
    df_info_result["GCA"] = df["gca_chain"].astype(str) + "." + df["gca_version"].astype(str)
    df_info_result["Reference genome"] = df_info_result["GCA"].apply(is_reference_genome)
    df_info_result["Release date"] = df["release_date"]

    # Load clade data
    clade_data = load_clade_data()

    # Add 'Clade' column to the info result by assigning clade based on lowest_taxon_id
    df_info_result["Internal clade"] = df_info_result["lowest_taxon_id"].apply(lambda x: assign_clade(x, clade_data))
    columns_to_drop_info = ['lowest_taxon_id']
    df_info_result.drop(columns=columns_to_drop_info, inplace=True, errors='ignore')


    # Clean genome_coverage by removing 'x' and converting to float
    df['metrics_value'] = df.apply(
        lambda row: float(row['metrics_value'].rstrip('x')) if row['metrics_name'] == 'genome_coverage' else row['metrics_value'], axis=1)

    # Pivot the data so each metric_name becomes a separate column and combine gca_chain and gca_version, correct date format
    df["GCA ID"] = df["gca_chain"].astype(str) + "." + df["gca_version"].astype(str)

    df_wide = df.pivot_table(index=["bioproject_id", "asm_level", "asm_type", "GCA ID", "release_date"], columns="metrics_name", values="metrics_value", aggfunc='first')

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

    # Apply `release_date` filter
    if release_date:
        release_date = pd.to_datetime(release_date)
        df_wide = df_wide[df_wide['release_date'] >= release_date]

    # Check if any assemblies meet the given thresholds
    if df_wide.empty:
        return "No assemblies meet the given thresholds.", refseq_count, None, None, None

    # Keep only rows in df_info_result that match the filtered GCA IDs
    filtered_gcas = df_wide["GCA ID"].unique()
    df_info_result = df_info_result[df_info_result["GCA"].isin(filtered_gcas)]

    # Drop specific columns
    columns_to_drop = ['contig_l50', 'gc_count', 'number_of_component_sequences', 'scaffold_l50', 'total_ungapped_length', 'number_of_organelles','total_number_of_chromosomes']  # Adjust this list as needed
    df_wide.drop(columns=columns_to_drop, inplace=True, errors='ignore')

    # Calculate summary statistics (Avg, Min, Max) for selected metrics
    summary_metrics = ['contig_n50', 'scaffold_n50', 'total_sequence_length', 'gc_percent', 'genome_coverage']
    summary_df = df_wide[summary_metrics].agg(['mean', 'min', 'max'])

    # Create the new GCA table
    df_gca_list = df_wide[["GCA ID"]]

    # Rename columns
    df_wide.rename(
        columns={'bioproject_id': 'BioProject ID', 'asm_level': 'Assembly level', 'number_of_contigs': 'Number of contigs', 'number_of_scaffolds': 'Number of scaffolds', 'scaffold_n50': 'Scaffold N50', 'total_sequence_length':'Sequence length', 'GCA ID': "GCA", 'contig_n50': 'Contig N50', 'gc_percent': 'GC%', 'genome_coverage': 'Genome coverage X', 'asm_type': "Assembly type"}, inplace=True)
    summary_df.rename(columns={'genome_coverage': 'Genome coverage X', 'contig_n50': 'Conting N50', 'scaffold_n50': 'Scaffold N50',  'total_sequence_length': "Sequence length", 'gc_percent': 'GC%'}, inplace=True)
    df_info_result.rename(
        columns={'bioproject_id': 'BioProject ID', 'scientific_name': 'Scientific name', 'common_name': "Common name", 'Reference genome': "Reference genome", 'group_name': "Custom group", 'release_date': "Release date"}, inplace=True)


    return df_wide, refseq_count, summary_df, df_info_result, df_gca_list

def main():
    parser = argparse.ArgumentParser(description="Fetch filtered assemblies for a given BioProject.")
    parser.add_argument('--bioproject_id', type=str, nargs='+', required=True, help="One or more BioProject IDs")
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

    args = parser.parse_args()

    metric_thresholds = {k: v[0] if isinstance(v, list) else v for k, v in vars(args).items() if v is not None and k not in ["bioproject_id", "output_dir", "asm_level", "asm_type", "release_date"]}

    all_metrics = ["gc_percent", "total_sequence_length", "contig_n50", "number_of_contigs", "number_of_scaffolds", "scaffold_n50", "genome_coverage"]

    df, refseq_count, summary_df, info_result, df_gca_list = get_filtered_assemblies(args.bioproject_id, metric_thresholds, all_metrics, args.asm_level, args.asm_type, args.release_date)


    # Check if 'df' is a string (error message)
    if isinstance(df, str):
        print(df)  # Print the error message
        return

    if df.empty:
        print("No assemblies meet the given thresholds.")
        return

    print("Assemblies that meet the given thresholds:")
    print(df)
    print(f"\nNumber of assemblies with a RefSeq accession: {refseq_count}")
    print("\nSummary statistics (Avg/Min/Max) for selected metrics:")
    print(summary_df)
    print("Assembly info:")
    print(info_result)
    print("GCA list:")
    print(df_gca_list)

    # Ensure the output directory exists
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Save the result as a CSV file
    output_filename = os.path.join(args.output_dir, f"{'_'.join(args.bioproject_id)}_filtered_assemblies.csv")
    df.to_csv(output_filename, index=False)

    # Save the summary statistics as a CSV file
    summary_filename = os.path.join(args.output_dir, f"{'_'.join(args.bioproject_id)}_filtered_assemblies_summary_statistics.csv")
    summary_df.to_csv(summary_filename, index=False)

    # Save the assembly info as a CSV file
    info_result_filename = os.path.join(args.output_dir, f"{'_'.join(args.bioproject_id)}_filtered_assemblies_info_result.csv")
    info_result.to_csv(info_result_filename, index=False)

    # Save the GCA list in a file
    gca_list_filename = os.path.join(args.output_dir, f"{'_'.join(args.bioproject_id)}_filtered_assemblies_gca_list.csv")
    df_gca_list.to_csv(gca_list_filename, index=False, header=False)

    print(f"\nThe results have been saved to {output_filename}")
    print(f"\nSummary statistics have been saved to {summary_filename}")
    print(f"\nAssembly info has been saved to {info_result_filename}")
    print(f"\nGCA list has been saved to {gca_list_filename}")


if __name__ == "__main__":
    main()
