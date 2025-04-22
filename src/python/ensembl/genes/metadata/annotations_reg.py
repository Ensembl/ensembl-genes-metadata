import os
import json
import time
import pymysql
import logging
import asyncio
import argparse
import requests
import pandas as pd
from datetime import datetime
from assemblies import get_filtered_assemblies
from assemblies import get_descendant_taxa
from check_transcriptomic_data import check_data_from_ena

# Configure logging
# Setup logging configuration
def setup_logging(output_dir):
    log_filename = os.path.join(output_dir, "annotations_log.txt")
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


def load_db_config():
    """Load database credentials from external JSON file."""
    config_path = "conf/db_config.json"
    if not os.path.exists(config_path):
        logging.error("Config file not found")
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


def load_bioproject_mapping():
    """Load database credentials from external JSON file."""
    bioproject_mapping = "data/bioproject_mapping.json"
    if not os.path.exists(bioproject_mapping):
        logging.error("BioProject mapping JSON not found")
        raise FileNotFoundError(f"BioProject mapping JSON '{bioproject_mapping}' not found.")

    with open(bioproject_mapping, "r") as f:
        return json.load(f)


def get_annotation(taxon_id, bioproject_id, annotation_release_date, annotation_date):
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
    parameters = []

    if bioproject_id:
        conditions.append(f"b.bioproject_id IN ({','.join(['%s'] * len(bioproject_id))})")
        parameters.extend(bioproject_id)
        logging.info(f"Filtering by BioProject IDs: {', '.join(bioproject_id)}")

    if taxon_id:
        descendant_taxa = get_descendant_taxa(taxon_id)
        logging.info(f"Retrieving annotation for taxon ID {taxon_id}.")
        if not descendant_taxa:
            return f"No descendant taxa found for Taxon ID {taxon_id}.", None

        conditions.append(f"a.lowest_taxon_id IN ({','.join(['%s'] * len(descendant_taxa))})")
        parameters.extend(descendant_taxa)

    if annotation_date:
        conditions.append("gb.date_started >= %s")
        parameters.append(annotation_date)
        logging.info(f"Filtering by release date: {', '.join(annotation_date)}")


    # If there are conditions, join them with AND; otherwise, select all
    where_clause = " WHERE " + " AND ".join(conditions) if conditions else ""

    query = f"""
                SELECT b.bioproject_id, a.lowest_taxon_id, 
                    gb.gca_accession, gb.gb_status, gb.genebuilder, gb.annotation_source, gb.date_completed, gb.date_completed_beta, gb.release_type,
                    gb.release_date_beta, gb.release_version_beta, gb.annotation_method, gb.date_started, s.scientific_name
                FROM bioproject b
                JOIN assembly a ON b.assembly_id = a.assembly_id
                JOIN genebuild gb ON a.assembly_id = gb.assembly_id
                JOIN species s ON a.lowest_taxon_id = s.lowest_taxon_id
                LEFT JOIN group_assembly g ON a.assembly_id = g.assembly_id
                {where_clause}
                ORDER BY gb.gca_accession;
            """

    cursor.execute(query, parameters)
    results = cursor.fetchall()
    conn.close()

    df_reg = pd.DataFrame(results)


    logging.debug(f"Total records fetched from registry: {len(df_reg)}")


    # Fill missing 'release_date_beta' with a placeholder value
    #df_reg['release_date_beta'] = df_reg['release_date_beta'].fillna('Not yet in Beta')


    # Apply release date filter only if provided
    if annotation_release_date:
	    annotation_release_date = pd.to_datetime(annotation_release_date)
	    df_reg['release_date_beta'] = pd.to_datetime(df_reg['release_date_beta'], errors='coerce')
	    df_reg = df_reg[df_reg['release_date_beta'] >= annotation_release_date]
	    logging.info(f"Showing records after {annotation_release_date}.")

    return df_reg


def summarize_assemblies(df):
    """Summarize the data by assembly level"""
    # Bin by assembly level
    level_summary = df.groupby('Assembly level',observed=False).size().reset_index(name='Number of Assemblies')

    return level_summary

def bin_by_genebuild_method(df):
    """Bins assemblies based on genebuild.method."""
    if 'annotation_method' not in df.columns:
        logging.error("Column 'annotation_method' not found in dataset.")
        return pd.DataFrame()  # Return empty DataFrame if column is missing

    # Group by genebuild.method and count occurrences
    method_summary = df.groupby('annotation_method', observed=False).size().reset_index(name='Number of Assemblies')

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

    # Extract annotation year from 'date_completed'
    if 'date_completed' in filtered_df.columns:
        filtered_df['year'] = pd.to_datetime(filtered_df['date_completed']).dt.year
        annotated_per_year = filtered_df.groupby('year').size().reset_index(name='Number of Annotated Genomes')
    else:
        logging.error("No date_completed column found in annotations.")
        return pd.DataFrame()
    logging.info("Created summary tables." )

    # Extract Ensembl release year
    if 'release_date_beta' in filtered_df.columns:
        filtered_df = filtered_df[filtered_df['release_date_beta'] != 'Not yet released']
        filtered_df = filtered_df.copy()  # Ensures we're working on a new copy
        filtered_df['beta_release_year'] = pd.to_datetime(filtered_df['release_date_beta'], errors='coerce').dt.year
        ensembl_release_summary = filtered_df.groupby('beta_release_year').size().reset_index(
            name='Number of Ensembl Releases')
    else:
        ensembl_release_summary = pd.DataFrame()

    # Merge all summaries
    yearly_summary = pd.merge(assemblies_per_year, annotated_per_year, on='year', how='outer').fillna(0)
    yearly_summary = pd.merge(yearly_summary, annotation_candidates, on='year', how='outer').fillna(0)
    yearly_summary = pd.merge(yearly_summary, ensembl_release_summary, left_on='year', right_on='beta_release_year', how='outer').fillna(0).drop(
        columns=['beta_release_year'])

    return yearly_summary


def bin_by_assembly_level(df):
    """Bins assemblies based on their assembly level."""
    if 'asm_level' not in df.columns:
        logging.error("Column 'Assembly level' not found in dataset.")
        return pd.DataFrame()  # Return empty DataFrame if column is missing

    # Group by Assembly level and count occurrences
    level_summary = df.groupby('asm_level', observed=False).size().reset_index(name='number_of_assemblies')
    logging.info("Binning complete." )
    return level_summary


def check_most_updated_annotation(df_info_result, filtered_df):
    """Checks if each annotated assembly is the latest available version and includes species name and latest GCA."""
    logging.info("Checking if annotations are the latest available GCA version." )
    # Extract version number from GCA identifier
    asm_info = df_info_result
    asm_info['version'] = asm_info['gca'].str.extract(r'GCA_\d+\.(\d+)').astype(float)
    asm_info['gca_latest'] = asm_info['gca']
    ann_filtered = filtered_df
    ann_filtered['version'] = ann_filtered['gca_accession'].str.extract(r'GCA_\d+\.(\d+)').astype(float)
    ann_filtered['gca_annotated'] = ann_filtered['gca_accession']

    # Create a DataFrame with all GCA and Version pairs
    latest_versions = asm_info[['gca', 'version', 'gca_latest', 'lowest_taxon_id', 'species_taxon_id']].rename(
        columns={'version': 'latest_version'}
    )
    latest_versions['gca'] = latest_versions['gca'].str.replace(r'\.\d+$', '', regex=True)

    # Rename GCA column in filtered_df for clarity
    ann_filtered['gca'] = ann_filtered['gca_annotated'].str.replace(r'\.\d+$', '', regex=True)

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

    merged_df = merged_df.rename(columns={'lowest_taxon_id_x': 'lowest_taxon_id'})

    return merged_df[['scientific_name', 'lowest_taxon_id', 'species_taxon_id', 'gca_annotated', 'gca_latest', 'latest_annotated']]


def main():
    parser = argparse.ArgumentParser(description='Fetches annotations and assemblies')
    parser.add_argument("--release_date", type=str, default="2000-01-01", help="Release date in YYYY-MM-DD format (default: 2000-01-01).")
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
    parser.add_argument('--current', type=int, choices=[0, 1], default=0, help="Check if GCA is the most current version (1 for yes, 0 for no)")
    parser.add_argument("--annotation_release_date", type=str, default="2000-01-01",
                        help="Annotation release date in YYYY-MM-DD format (default: 2000-01-01).")
    parser.add_argument("--annotation_date", type=str, default="2000-01-01", help="Annotation date in YYYY-MM-DD format (default: 2000-01-01).")

    args = parser.parse_args()

    release_date = args.release_date
    output_dir = args.output_dir
    taxon_id = args.taxon_id
    asm_level = args.asm_level
    asm_type = args.asm_type
    bioproject_id = args.bioproject_id
    annotation_release_date = args.annotation_release_date
    annotation_date = args.annotation_date

    # Ensure the output directory exists
    if not os.path.exists(args.output_dir):
        logging.info(f"Creating output directory {args.output_dir}")
        os.makedirs(args.output_dir)

    # Set up logging with the user-specified output directory
    setup_logging(args.output_dir)

    # Define metrics to consider
    metric_thresholds = {k: v[0] if isinstance(v, list) else v for k, v in vars(args).items() if
                         v is not None and k not in ["bioproject_id", "output_dir", "asm_level", "asm_type",
                                                     "release_date", "reference", "taxon_id","current"]}

    all_metrics = ["gc_percent", "total_sequence_length", "contig_n50", "number_of_contigs", "number_of_scaffolds",
                   "scaffold_n50", "genome_coverage"]

    # Fetch assemblies
    df_wide, summary_df, df_info_result, df_gca_list, taxonomy_dict = get_filtered_assemblies(bioproject_id, metric_thresholds, all_metrics, asm_level, asm_type, release_date, taxon_id, args.current)
    logging.info(f"Filtered assemblies: {len(df_wide)}")


    # Check for transcriptomic data only if the user requested it
    if args.trans:
        # Check transcriptomic data for each taxon_id in the dataset
        logging.info(f"Transcriptomic data check requested requested")
        # Get unique taxon IDs and filter out NaN values
        taxon_ids = [tid for tid in df_info_result["lowest_taxon_id"].unique() if pd.notna(tid)]
        species_taxon_ids = [tid for tid in df_info_result["species_taxon_id"].unique() if pd.notna(tid)]
        genus_taxon_ids = [gtid for gtid in df_info_result["genus_taxon_id"].unique() if pd.notna(gtid)]

        # Count NaN values and log warning if any are found
        nan_lowest_count = df_info_result["lowest_taxon_id"].isna().sum()
        nan_species_count = df_info_result["species_taxon_id"].isna().sum()
        nan_genus_count = df_info_result["genus_taxon_id"].isna().sum()

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
        df_info_result = df_info_result.merge(
            transcriptomic_df, left_on="lowest_taxon_id", right_on="Taxon ID", how="left", suffixes=('', '_lowest')
        )

        # Merge for the species taxon ID
        df_info_result = df_info_result.merge(
            transcriptomic_df, left_on="species_taxon_id", right_on="Taxon ID", how="left",
            suffixes=('_lowest', '_species')
        )

        # Merge for the genus taxon ID (separate column)
        df_info_result = df_info_result.merge(
            transcriptomic_df, left_on="genus_taxon_id", right_on="Taxon ID", how="left", suffixes=('_lowest', '_genus')
        )

        # Drop redundant 'Taxon ID' columns (both for lowest and genus)
        df_info_result.drop(columns=["Taxon ID_lowest", "Taxon ID_species", "Taxon ID"], inplace=True)


    if isinstance(df_wide, str):  # Check if there's an error message
        print(df_wide)
        logging.error(f"Error message: {df_wide})") # Print the error message
        return

    # Bin assemblies by Assembly level
    level_summary = bin_by_assembly_level(df_wide)

    # Fetch dataset annotations using the same release date
    filtered_df = get_annotation(taxon_id, bioproject_id, annotation_release_date, annotation_date)

    # Create the FTP URL using the scientific_name, replacing spaces with underscores
    logging.info("Generating FTP paths.")
    filtered_df['ftp'] = filtered_df.apply(
        lambda
            row: f"https://ftp.ebi.ac.uk/pub/ensemblorganisms/{row['scientific_name'].replace(' ', '_')}/{row['gca_accession']}/"
        if pd.notnull(row['scientific_name']) and pd.notnull(row['gca_accession']) and pd.notnull(
            row.get('beta_release_date')) else None,
        axis=1)

    # Bin the filtered results by genebuild.method
    method_summary = bin_by_genebuild_method(filtered_df)

    yearly_summary = generate_yearly_summary(df_wide, df_info_result, filtered_df)

    #Check annotation statusin beta and main
    logging.info("Checking assemblies in annotation table.")
    df_info_result['annotated'] = df_info_result['gca'].isin(
        filtered_df['gca_accession']).map({True: 'Yes', False: 'No'})

    #Check if annotated is the latest GCA version in the registry
    latest_annotated_df = check_most_updated_annotation(df_info_result, filtered_df)

    #filtered_df = filtered_df.drop(columns=['year', 'gca', 'version'])
    #df_info_result = df_info_result.drop(columns=['year', 'version', 'gca_latest'])

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
    logging.info("Finished saving tables.")

if __name__ == "__main__":
    main()
