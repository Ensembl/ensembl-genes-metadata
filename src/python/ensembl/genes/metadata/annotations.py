import os
import json
import time
import pymysql
import logging
import asyncio
import argparse
import requests
import pandas as pd
from typing import List
from datetime import datetime
from assemblies import get_filtered_assemblies
from check_transcriptomic_data import check_data_from_ena

# Configure logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def load_db_config():
    """Load database credentials from external JSON file."""
    config_path = "conf/prod_dbs_conf.json"
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Database config file '{config_path}' not found.")

    with open(config_path, "r") as f:
        return json.load(f)


def connect_db(config_key):
    """Establish connection to the MySQL database."""
    db_config = load_db_config()
    return pymysql.connect(**db_config[config_key], cursorclass=pymysql.cursors.DictCursor)


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
        "retmax": 100000,  # Fetch in chunks this is the max value
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

def get_gca_accessions(bioproject_id) -> List[str]:
    """
    Fetches GCA assembly accessions from a given NCBI BioProject ID using the NCBI Datasets API.

    Args:
        bioproject_id (str): The NCBI BioProject ID (e.g., "PRJEB40665")

    Returns:
        list: A list of GCA accession numbers associated with the BioProject
    """
    base_url = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/bioproject"
    next_page_token = None
    gca_accessions = []

    while True:
        url = f"{base_url}/{bioproject_id}/dataset_report"
        params = {}
        if next_page_token:
            params['page_token'] = next_page_token

        try:
            response = requests.get(url, params=params)
            response.raise_for_status()  # This will raise an exception for 4XX and 5XX errors
            data = response.json()

            assemblies = data.get('reports', [])
            for assembly in assemblies:
                assembly_accession = assembly.get('accession')

                # Only include GCA accessions (GenBank)
                if assembly_accession and assembly_accession.startswith('GCA_'):
                    gca_accessions.append(assembly_accession)

            # Check for the next page token
            next_page_token = data.get('next_page_token')
            if not next_page_token:
                break

            # Inside the loop where you make requests
            if response.status_code == 429:
                print("Rate limit exceeded, sleeping for 10 seconds...")
                time.sleep(10)
                continue  # retry the request after delay

        except requests.HTTPError as http_err:
            print(f"HTTP error occurred: {http_err}")
            break
        except Exception as err:
            print(f"An error occurred: {err}")
            break

    logging.debug(f"Final GCA accessions: {gca_accessions}")
    return gca_accessions


def get_annotation(release_date, taxon_id, bioproject_id):
    db_config_key = "prod"
    conn = connect_db(db_config_key)
    cursor = conn.cursor()

    if isinstance(bioproject_id, list) and len(bioproject_id) > 0:
        bioproject_id = bioproject_id[0]  # Take the first Bioproject ID

    print(f"Bioproject ID: {bioproject_id}")

    # Fetch GCA list from NCBI if BioProject ID is provided
    gca_accessions = get_gca_accessions(bioproject_id) if bioproject_id else None
    logging.debug(f"GCA Accessions: {gca_accessions}")


    # Build dynamic SQL filtering
    conditions = []
    parameters = []

    if taxon_id:
        descendant_taxa = get_descendant_taxa(taxon_id)
        if not descendant_taxa:
            return f"No descendant taxa found for Taxon ID {taxon_id}.", None

        conditions.append(f"o.species_taxonomy_id IN ({','.join(['%s'] * len(descendant_taxa))})")
        parameters.extend(descendant_taxa)

    if gca_accessions:  # If filtering by BioProject ID, add GCA condition
        conditions.append(f"a.accession IN ({','.join(['%s'] * len(gca_accessions))})")
        parameters.extend(gca_accessions)

    # If there are conditions, join them with AND; otherwise, select all
    where_clause = " WHERE " + " AND ".join(conditions) if conditions else ""
    logging.debug(f"Params passed to the query: {parameters}")

    query = f"""
                SELECT da.value, da.attribute_id, a.accession AS gca, r.release_date AS ensembl_release_date, o.scientific_name, o.species_taxonomy_id
                FROM dataset_attribute da
                JOIN dataset d ON da.dataset_id = d.dataset_id
                JOIN genome_dataset gd ON d.dataset_id = gd.dataset_id
                JOIN genome g ON gd.genome_id = g.genome_id
                JOIN assembly a ON g.assembly_id = a.assembly_id
                JOIN organism o ON g.organism_id = o.organism_id
                LEFT JOIN ensembl_release r ON r.release_id = gd.release_id
                {where_clause}
                AND d.name = 'genebuild'
                AND da.attribute_id IN (34, 37, 25, 183, 31, 40, 42, 48, 170, 56, 212);
            """

    cursor.execute(query, parameters)
    print(f"Params passed to the query: {parameters}")
    print(f"Final WHERE clause: {where_clause}")
    print(f"Executing query: {query}")
    results = cursor.fetchall()
    conn.close()

    df_prod = pd.DataFrame(results)

    logging.debug(f"Total records fetched from production db: {len(df_prod)}")

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
        48: "nc_small_non_coding_genes",
        170: "nc_total_exons",
        56: "ps_average_sequence_length",
    })

    # Fill missing 'ensembl_release_date' with a placeholder value
    df_prod['ensembl_release_date'] = df_prod['ensembl_release_date'].fillna('Not yet released')

    # Pivot to wide format without losing records
    df_pivoted = df_prod.pivot_table(
        index=["ensembl_release_date", "gca", "scientific_name", "species_taxonomy_id"],
        columns="attribute_id",
        values="value",
        aggfunc="first"  # or use 'list' if you want to capture multiple values for the same attribute_id
    ).reset_index()

    logging.debug(f"Pivoted DataFrame shape before filtering: {df_pivoted.shape}")


    # Ensure genebuild.last_geneset_update has the correct format by adding '-01' for day precision
    df_pivoted['last_geneset_update'] = pd.to_datetime(
        df_pivoted['last_geneset_update'] + '-01', errors='coerce')

    # Convert the release_date to a datetime object
    release_date = pd.to_datetime(release_date)

    # Apply the release_date filter to 'genebuild.last_geneset_update'
    df_filtered = df_pivoted[df_pivoted['last_geneset_update'] >= release_date]
    logging.debug(f"Records after release date filtering: {len(df_filtered)}")

    df_filtered = df_filtered[df_filtered['gca'].str.startswith('GCA')]
    logging.debug(f"Records after GCA filtering: {len(df_filtered)}")

    logging.info(f"Filtered {len(df_filtered)} records based on the release date.")

    return df_filtered


def summarize_assemblies(df):
    """Summarize the data by assembly level"""
    # Bin by assembly level
    level_summary = df.groupby('Assembly level',observed=False).size().reset_index(name='Number of Assemblies')

    return level_summary

def bin_by_genebuild_method(df):
    """Bins assemblies based on genebuild.method."""
    if 'genebuild_method' not in df.columns:
        print("Column 'genebuild.method' not found in dataset.")
        return pd.DataFrame()  # Return empty DataFrame if column is missing

    # Group by genebuild.method and count occurrences
    method_summary = df.groupby('genebuild_method', observed=False).size().reset_index(name='Number of Assemblies')

    return method_summary


def generate_yearly_summary(df_wide, df_info_result, filtered_df):
    """Generate a table summarizing assemblies above contig level and annotated genomes per year."""

    assembly_info = df_info_result
    assemblies = df_wide
    assembly_info['year'] = pd.to_datetime(assembly_info['release_date']).dt.year
    assembly_info = assemblies.merge(assembly_info[['gca', 'year']], on='gca', how='left')

    valid_levels = ["Scaffold", "Chromosome", "Complete genome"]

    # Filter only the valid assembly levels
    df_filtered_levels = assembly_info[assembly_info['asm_level'].isin(valid_levels)].copy()

    # Count assemblies per year
    assemblies_per_year = df_filtered_levels.groupby('year').size().reset_index(
        name='Number of Assemblies')

    # Ensure 'Contig N50' is numeric
    df_filtered_levels['contig_n50'] = pd.to_numeric(df_filtered_levels['contig_n50'], errors='coerce')


    # Filter for assemblies with contigN50 >= 100000 (Annotation Candidates)
    annotation_candidates = df_filtered_levels[df_filtered_levels['contig_n50'] >= 100000].groupby(
        'year').size().reset_index(
        name='Annotation Candidates')

    # Extract annotation year from 'last_geneset_update'
    if 'last_geneset_update' in filtered_df.columns:
        filtered_df['year'] = pd.to_datetime(filtered_df['last_geneset_update']).dt.year
        annotated_per_year = filtered_df.groupby('year').size().reset_index(name='Number of Annotated Genomes')
    else:
        print("No 'last_geneset_update' column found in annotations.")
        return pd.DataFrame()

    # Extract Ensembl release year
    if 'ensembl_release_date' in filtered_df.columns:
        filtered_df = filtered_df[filtered_df['ensembl_release_date'] != 'Not yet released']
        filtered_df = filtered_df.copy()  # Ensures we're working on a new copy
        filtered_df['ensembl_release_year'] = pd.to_datetime(filtered_df['ensembl_release_date'], errors='coerce').dt.year
        ensembl_release_summary = filtered_df.groupby('ensembl_release_year').size().reset_index(
            name='Number of Ensembl Releases')
    else:
        ensembl_release_summary = pd.DataFrame()

    # Merge all summaries
    yearly_summary = pd.merge(assemblies_per_year, annotated_per_year, on='year', how='outer').fillna(0)
    yearly_summary = pd.merge(yearly_summary, annotation_candidates, on='year', how='outer').fillna(0)
    yearly_summary = pd.merge(yearly_summary, ensembl_release_summary, left_on='year', right_on='ensembl_release_year', how='outer').fillna(0).drop(
        columns=['ensembl_release_year'])

    return yearly_summary


def bin_by_assembly_level(df):
    """Bins assemblies based on their assembly level."""
    if 'asm_level' not in df.columns:
        print("Column 'Assembly level' not found in dataset.")
        return pd.DataFrame()  # Return empty DataFrame if column is missing

    # Group by Assembly level and count occurrences
    level_summary = df.groupby('asm_level', observed=False).size().reset_index(name='number_of_assemblies')

    return level_summary


def check_most_updated_annotation(df_info_result, filtered_df):
    """Checks if each annotated assembly is the latest available version and includes species name and latest GCA."""

    # Extract version number from GCA identifier
    asm_info = df_info_result
    asm_info['version'] = asm_info['gca'].str.extract(r'GCA_\d+\.(\d+)').astype(float)
    asm_info['gca_latest'] = asm_info['gca']
    ann_filtered = filtered_df
    ann_filtered['version'] = ann_filtered['gca'].str.extract(r'GCA_\d+\.(\d+)').astype(float)
    ann_filtered['gca_annotated'] = ann_filtered['gca']

    # Create a DataFrame with all GCA and Version pairs
    latest_versions = asm_info[['gca', 'version', 'gca_latest', 'lowest_taxon_id', 'species_taxon_id']].rename(
        columns={'version': 'latest_version'}
    )
    latest_versions['gca'] = latest_versions['gca'].str.replace(r'\.\d+$', '', regex=True)

    # Rename GCA column in filtered_df for clarity
    ann_filtered['gca'] = ann_filtered['gca'].str.replace(r'\.\d+$', '', regex=True)

    # Merge latest versions with filtered_df based on GCA
    merged_df = ann_filtered.merge(latest_versions, on="gca", how='left')

    # Create new columns for annotated and assembly versions
    merged_df['annotated_version'] = merged_df['version']
    merged_df['assembly_version'] = merged_df['latest_version']


    # Check if the annotated GCA is the same as the latest assembly GCA
    def check_latest_annotated(row):
        if pd.isna(row['assembly_version']):
            return 'Yes, low quality assembly version'
        return 'Yes' if row['annotated_version'] == row['assembly_version'] else 'No'

    merged_df['latest_annotated'] = merged_df.apply(check_latest_annotated, axis=1)

    # Order the DataFrame by Scientific name
    merged_df.sort_values(by='scientific_name', inplace=True)
    print("merged df for debug:")
    print(merged_df)

    return merged_df[['scientific_name', 'lowest_taxon_id', 'species_taxon_id', 'gca_annotated', 'gca_latest', 'latest_annotated']]



def main():
    parser = argparse.ArgumentParser(description='Fetches annotations and assemblies')
    parser.add_argument("--release_date", type=str, default="2000-01-01", help="Release date in YYYY-MM-DD format (default: 2019-01-01).")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save output files.")
    parser.add_argument("--taxon_id", type=str, required=False, help="Taxon ID for filtering (e.g., 40674 for Mammalia).")
    parser.add_argument('--asm_level', type=str, nargs='+', help="Assembly level options: 'Contig', 'Scaffold', 'Chromosome', 'Complete genome'.")
    parser.add_argument('--asm_type', type=str, nargs='+', help="Assembly type: 'haploid', 'alternate-pseudohaplotype', 'unresolved-diploid', 'haploid-with-alt-loci', 'diploid'.")
    parser.add_argument('--contig_n50', type=float, help="Contig N50")
    parser.add_argument('--gc_percent', type=float, help="GC percent threshold")
    parser.add_argument('--total_sequence_length', type=float, help="Total sequence length in bp")
    parser.add_argument('--number_of_contigs', type=float, help="Number of contigs")
    parser.add_argument('--number_of_scaffolds', type=float, help="Number of scaffolds")
    parser.add_argument('--scaffold_n50', type=float, help="Scaffold N50")
    parser.add_argument('--genome_coverage', type=float, help="Genome coverage in bp")
    parser.add_argument('--bioproject_id', type=str, nargs='+', help="One or more BioProject IDs")
    parser.add_argument('--trans', type=int, choices=[0, 1], default=0, help="Check if a taxon ID has transcriptomic data from ENA (1 for yes, 0 for no)")

    args = parser.parse_args()

    release_date = args.release_date
    output_dir = args.output_dir
    taxon_id = args.taxon_id
    asm_level = args.asm_level
    asm_type = args.asm_type
    bioproject_id = args.bioproject_id

    os.makedirs(output_dir, exist_ok=True)

    # Define metrics to consider
    metric_thresholds = {k: v[0] if isinstance(v, list) else v for k, v in vars(args).items() if
                         v is not None and k not in ["bioproject_id", "output_dir", "asm_level", "asm_type",
                                                     "release_date", "reference", "taxon_id"]}

    all_metrics = ["gc_percent", "total_sequence_length", "contig_n50", "number_of_contigs", "number_of_scaffolds",
                   "scaffold_n50", "genome_coverage"]

    # Fetch assemblies released in the past 5 years
    df_wide, summary_df, df_info_result, df_gca_list = get_filtered_assemblies(bioproject_id, metric_thresholds, all_metrics, asm_level, asm_type, release_date, taxon_id)


    # Check for transcriptomic data only if the user requested it
    if args.trans:
        # Check transcriptomic data for each taxon_id in the dataset
        taxon_ids = df_info_result["lowest_taxon_id"].unique()
        genus_taxon_ids = df_info_result["genus_taxon_id"].unique()
        all_taxon_ids = set(taxon_ids).union(set(genus_taxon_ids))

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
        df_info_result = df_info_result.merge(
            transcriptomic_df, left_on="lowest_taxon_id", right_on="Taxon ID", how="left"
        )

        # Merge for the genus taxon ID (separate column)
        df_info_result = df_info_result.merge(
            transcriptomic_df, left_on="genus_taxon_id", right_on="Taxon ID", how="left", suffixes=('_lowest', '_genus')
        )

        # Drop redundant 'Taxon ID' columns (both for lowest and genus)
        df_info_result.drop(columns=["Taxon ID_lowest", "Taxon ID_genus"], inplace=True)

        # Rename or reorder columns if needed, to better organize your data
        #df_info_result.rename(columns={"Taxon ID": "Transcriptional Data"}, inplace=True)

    if isinstance(df_wide, str):  # Check if there's an error message
        print(df_wide)
        return

    # Bin assemblies by Assembly level
    level_summary = bin_by_assembly_level(df_wide)

    # Fetch dataset annotations using the same release date
    filtered_df = get_annotation(release_date, taxon_id, bioproject_id)

    # Create the FTP URL using the scientific_name, replacing spaces with underscores
    filtered_df['ftp'] = filtered_df.apply(
        lambda
            row: f"https://ftp.ebi.ac.uk/pub/ensemblorganisms/{row['scientific_name'].replace(' ', '_')}/{row['gca']}/"
        if pd.notnull(row['scientific_name']) and pd.notnull(row['gca']) and pd.notnull(
            row.get('ensembl_release_date')) else None,
        axis=1)

    # Bin the filtered results by genebuild.method
    method_summary = bin_by_genebuild_method(filtered_df)

    yearly_summary = generate_yearly_summary(df_wide, df_info_result, filtered_df)

    #Check annotation status
    # Add an 'Annotated' column to df_info_result
    df_info_result['annotated'] = df_info_result['gca'].str.extract(r'(GCA_\d+)', expand=False).isin(
        filtered_df['gca'].str.extract(r'(GCA_\d+)', expand=False)).map({True: 'Yes', False: 'No'})

    #Check if annotated is the latest GCA version in the registry
    latest_annotated_df = check_most_updated_annotation(df_info_result, filtered_df)

    filtered_df = filtered_df.drop(columns=['year', 'gca', 'version'])
    df_info_result = df_info_result.drop(columns=['year', 'version', 'gca_latest'])

    print("\nDataset filtered:")
    print(method_summary)

    print("\nDataset yearly summary:")
    print(yearly_summary)

    print("\nAnnotated:")
    print(filtered_df)

    print("\nAssembly Level Summary:")
    print(level_summary)

    print("\nMost Updated Assemblies Annotated:")
    print(latest_annotated_df)

# Get current date in YYYYMMDD format
    date_str = datetime.now().strftime("%Y%m%d")

# Save results to CSV files
    method_summary.to_csv(os.path.join(output_dir, f"annotation_method_summary_{date_str}.csv"), index=False)
    yearly_summary.to_csv(os.path.join(output_dir, f"yearly_summary_{date_str}.csv"), index=False)
    filtered_df.to_csv(os.path.join(output_dir, f"annotations_{date_str}.csv"), index=False)
    df_wide.to_csv(os.path.join(output_dir, f"assemblies_{date_str}.csv"), index=False)
    df_info_result.to_csv(os.path.join(output_dir, f"assembly_annotation_status_{date_str}.csv"), index=False)
    level_summary.to_csv(os.path.join(output_dir, f"assembly_level_summary_{date_str}.csv"), index=False)
    latest_annotated_df.to_csv(os.path.join(output_dir, f"annotation_GCA_update_status_{date_str}.csv"), index=False)

if __name__ == "__main__":
    main()
