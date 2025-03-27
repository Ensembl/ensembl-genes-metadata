import pandas as pd
import json
import pymysql
import argparse
import os
import requests
import asyncio
import time
from check_transcriptomic_data import check_data_from_ena

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


def get_taxonomy_from_db(lowest_taxon_id):
    """Retrieve taxonomy hierarchy from the MySQL database."""
    conn = connect_db()
    cursor = conn.cursor()

    query = """
    SELECT taxon_class_id, taxon_class 
    FROM taxonomy
    WHERE lowest_taxon_id = %s
    ORDER BY FIELD(taxon_class, 'species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom');
    """

    cursor.execute(query, (lowest_taxon_id,))
    taxonomy_hierarchy = cursor.fetchall()

    if not taxonomy_hierarchy:
        return None

    return taxonomy_hierarchy


def assign_clade_and_species(lowest_taxon_id, clade_data, chordata_taxon_id=7711):
    """Assign internal clade and species taxon ID based on taxonomy using the provided clade data,
       and check if the taxon ID is a descendant of the chordata taxon ID (7711)."""
    # Retrieve the taxonomy hierarchy for the given lowest_taxon_id
    taxonomy_hierarchy = get_taxonomy_from_db(lowest_taxon_id)

    if not taxonomy_hierarchy:
        return "Taxon hierarchy not found.", None

    species_taxon_id = None
    genus_taxon_id = None
    internal_clade = "Unassigned"  # Default value if no clade is found

    # Define the taxon classes in hierarchical order
    taxon_classes_order = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']

    # Loop through the taxonomy hierarchy in the order from species to kingdom
    for taxon_class in taxon_classes_order:
        # Find the taxon of the current class in the hierarchy
        matching_taxon = next((t for t in taxonomy_hierarchy if t['taxon_class'] == taxon_class), None)

        if matching_taxon:
            taxon_class_id = matching_taxon['taxon_class_id']  # Extract taxon_class_id

            # If it's the species level, always assign species_taxon_id
            if taxon_class == 'species' and species_taxon_id is None:
                species_taxon_id = taxon_class_id

            # If it's the genus level, always assign genus_taxon_id
            if taxon_class == 'genus' and genus_taxon_id is None:
                genus_taxon_id = taxon_class_id

            # Check for matching taxon_id in clade settings
            for clade_name, details in clade_data.items():
                if details.get("taxon_id") == taxon_class_id:
                    internal_clade = clade_name
                    break  # If a clade is found, we stop looking

            # If a clade is found, we stop looking
            if internal_clade != "Unassigned":
                break
    # Check if the taxon ID is 7742 or a descendant of it
    # Loop through the hierarchy and check if chordata_taxon_id (7711) is part of it
    if any(t['taxon_class_id'] == chordata_taxon_id for t in taxonomy_hierarchy):
        pipeline = "main"  # Assign "main" if it is 7742 or a descendant
    else:
        pipeline = "anno"  # Otherwise, assign "anno"

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
            print(f"Error retrieving taxonomic data from NCBI. HTTP {response.status_code}")
            break

        try:
            result = response.json()
            batch_ids = result.get("esearchresult", {}).get("idlist", [])
            if not batch_ids:
                break  # No more results

            taxon_ids.update(batch_ids)

            # Update retstart to fetch the next batch
            params["retstart"] += params["retmax"]

            # Respect NCBI rate limits
            time.sleep(0.5)  # Avoid overloading NCBI servers

        except Exception as e:
            print(f"Error processing response: {e}")
            break

    return taxon_ids

def get_filtered_assemblies(bioproject_id, metric_thresholds, all_metrics, asm_level, asm_type, release_date,taxon_id):
    """Fetch all assemblies and their metrics filter results based on given thresholds,
    and format the results with metrics as separate columns."""

    conn = connect_db()
    cursor = conn.cursor()

    # Check if bioproject_id is provided before validating
    if bioproject_id:
        cursor.execute("SELECT DISTINCT bioproject_id FROM bioproject;")
        valid_bioprojects = {row['bioproject_id'] for row in cursor.fetchall()}
        invalid_bioprojects = set(bioproject_id) - valid_bioprojects
        if invalid_bioprojects:
            conn.close()
            return f"The following BioProject IDs were not found: {', '.join(invalid_bioprojects)}", None, None, None

    # Build dynamic SQL filtering
    conditions = []
    params = []

    if bioproject_id:
        conditions.append(f"b.bioproject_id IN ({','.join(['%s'] * len(bioproject_id))})")
        params.extend(bioproject_id)

    if release_date:
        conditions.append("a.release_date >= %s")
        params.extend(release_date)

    if taxon_id:
        descendant_taxa = get_descendant_taxa(taxon_id)
        if not descendant_taxa:
            return f"No descendant taxa found for Taxon ID {taxon_id}.", None

        conditions.append(f"s.lowest_taxon_id IN ({','.join(['%s'] * len(descendant_taxa))})")
        params.extend(descendant_taxa)


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
    conn.close()

    # Convert results to a Pandas DataFrame
    df = pd.DataFrame(results)
    if df.empty:
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


    # Apply assembly level an assembly type filtering (supporting multiple values)
    if asm_level:
        df_wide = df_wide[df_wide['asm_level'].isin(asm_level)]

    if asm_type:
        df_wide = df_wide[df_wide['asm_type'].isin(asm_type)]

    # Check if any assemblies meet the given thresholds
    if df_wide.empty:
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
    df_info_result[['internal_clade', 'species_taxon_id', 'genus_taxon_id', 'pipeline']] = df_info_result['lowest_taxon_id'].apply(lambda x: pd.Series(assign_clade_and_species(x, clade_data)))

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
    parser.add_argument('--asm_level', type=str, nargs='+', help="Assembly level options: 'Contig' 'Scaffold' 'Chromosome' 'Complete genome'.")
    parser.add_argument('--asm_type', type=str, nargs='+', help="Assembly type: 'haploid' 'alternate-pseudohaplotype' 'unresolved-diploid' 'haploid-with-alt-loci' 'diploid'.")
    parser.add_argument('--output_dir', type=str, default='./', help="Directory to save the CSV file. Files will only be saved if assemblies meet the given thresholds.")
    parser.add_argument('--taxon_id', type=int, help="NCBI Taxon ID to filter by (e.g., 40674 for Mammalia)")
    parser.add_argument('--reference', type=int, choices=[0, 1], default=0, help="Check if GCA is a reference genome (1 for yes, 0 for no)")
    parser.add_argument('--trans', type=int, choices=[0, 1], default=0, help="Check if a taxon ID has transcriptomic data from ENA (1 for yes, 0 for no)")

    args = parser.parse_args()

    metric_thresholds = {k: v[0] if isinstance(v, list) else v for k, v in vars(args).items() if v is not None and k not in ["bioproject_id", "output_dir", "asm_level", "asm_type", "release_date", "reference", "taxon_id"]}

    all_metrics = ["gc_percent", "total_sequence_length", "contig_n50", "number_of_contigs", "number_of_scaffolds", "scaffold_n50", "genome_coverage"]

    df, summary_df, info_result, df_gca_list = get_filtered_assemblies(args.bioproject_id, metric_thresholds, all_metrics, args.asm_level, args.asm_type, args.release_date, args.taxon_id)

    # Check for reference genome only if the user requested it
    if args.reference:
        info_result["reference_genome"] = info_result["gca"].apply(is_reference_genome)

    # Check for transcriptomic data only if the user requested it
    if args.trans:
        # Check transcriptomic data for each taxon_id in the dataset
        taxon_ids = info_result["lowest_taxon_id"].unique()
        semaphore = asyncio.Semaphore(5)

        # Run transcriptomic data retrieval asynchronously
        async def fetch_transcriptomic_data():
            return await asyncio.gather(*[check_data_from_ena(taxon_id, tree=True, semaphore=semaphore) for taxon_id in taxon_ids])

        transcriptomic_results = asyncio.run(fetch_transcriptomic_data())

        # Convert results into a DataFrame
        transcriptomic_df = pd.DataFrame(transcriptomic_results)

        # Merge the transcriptomic data into info_result using 'Taxon ID' as the key
        info_result = info_result.merge(transcriptomic_df, left_on="lowest_taxon_id", right_on="Taxon ID", how="left")

        # Drop redundant 'Taxon ID' column
        info_result.drop(columns=["Taxon ID"], inplace=True)

    # Check if 'df' is a string (error message)
    if isinstance(df, str):
        print(df)  # Print the error message
        return

    if df.empty:
        print("No assemblies meet the given thresholds.")
        return

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
    bioproject_name = '_'.join(args.bioproject_id) if args.bioproject_id else ''
    release_date_name = args.release_date.replace('-', '') if args.release_date else ''

    # Create the output filename based on the provided values
    output_filename = os.path.join(args.output_dir,
                                   f"{bioproject_name}{'_' if bioproject_name and release_date_name else ''}{release_date_name}_filtered_assemblies.csv")
    df.to_csv(output_filename, index=False)

    # Save the summary statistics as a CSV file
    summary_filename = os.path.join(args.output_dir,
                                   f"{bioproject_name}{'_' if bioproject_name and release_date_name else ''}{release_date_name}_filtered_assemblies_summary_statistics.csv")
    summary_df.to_csv(summary_filename, index=False)

    # Save the assembly info as a CSV file
    info_result_filename = os.path.join(args.output_dir,
                                   f"{bioproject_name}{'_' if bioproject_name and release_date_name else ''}{release_date_name}_filtered_assemblies_info_result.csv")
    info_result.to_csv(info_result_filename, index=False)

    # Save the GCA list in a file
    gca_list_filename = os.path.join(args.output_dir,
                                        f"{bioproject_name}{'_' if bioproject_name and release_date_name else ''}{release_date_name}_filtered_assemblies_gca_list.csv")
    df_gca_list.to_csv(gca_list_filename, index=False, header=False)

    print(f"\nThe results have been saved to {output_filename}")
    print(f"\nSummary statistics have been saved to {summary_filename}")
    print(f"\nAssembly info has been saved to {info_result_filename}")
    print(f"\nGCA list has been saved to {gca_list_filename}")


if __name__ == "__main__":
    main()
