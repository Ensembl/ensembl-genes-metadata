import pandas as pd
import json
import pymysql
import os
import argparse
from GB_metadata_reporting import get_filtered_assemblies

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


def get_annotation(release_date):
    db_config_key = "prod"
    conn = connect_db(db_config_key)
    cursor = conn.cursor()

    query = f"""
                SELECT da.value, da.attribute_id, d.label, r.release_date AS ensembl_release_date
                FROM dataset_attribute da
                JOIN dataset d ON da.dataset_id = d.dataset_id
                JOIN genome_dataset g ON d.dataset_id = g.dataset_id
                LEFT JOIN ensembl_release r ON r.release_id = g.release_id
                WHERE da.attribute_id IN (34, 37, 25, 183, 31, 40, 42, 48, 170, 56, 212)
                AND d.name = 'genebuild'
            """

    cursor.execute(query)
    results = cursor.fetchall()
    conn.close()

    df_prod = pd.DataFrame(results)

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


    # Pivot to wide format without losing records
    df_pivoted = df_prod.pivot(index=["label", "ensembl_release_date"], columns="attribute_id", values="value").reset_index()
    # Extract only the GCA accession
    df_pivoted["GCA"] = df_pivoted["label"].str.extract(r"(GCA_\d+\.\d+)")
    df_pivoted = df_pivoted[df_pivoted['GCA'].notna() & (df_pivoted['GCA'] != '')]
    df_pivoted = df_pivoted.drop(columns=['label'])

    # Ensure genebuild.last_geneset_update has the correct format by adding '-01' for day precision
    df_pivoted['last_geneset_update'] = pd.to_datetime(
        df_pivoted['last_geneset_update'] + '-01', errors='coerce')

    # Convert the release_date to a datetime object
    release_date = pd.to_datetime(release_date)

    # Apply the release_date filter to 'genebuild.last_geneset_update'
    df_filtered = df_pivoted[df_pivoted['last_geneset_update'] >= release_date]

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


def generate_yearly_summary(df, df_info_result, filtered_df):
    """Generate a table summarizing assemblies above contig level and annotated genomes per year."""
    df_info_result['Year'] = pd.to_datetime(df_info_result['Release date']).dt.year
    df = df.merge(df_info_result[['GCA', 'Year']], on='GCA', how='left')

    valid_levels = ["Scaffold", "Chromosome", "Complete genome"]

    # Filter only the valid assembly levels
    df_filtered_levels = df[df['Assembly level'].isin(valid_levels)]

    # Count assemblies per year
    assemblies_per_year = df_filtered_levels.groupby('Year').size().reset_index(
        name='Number of Assemblies')

    # Ensure 'Contig N50' is numeric
    df_filtered_levels['Contig N50'] = pd.to_numeric(df_filtered_levels['Contig N50'], errors='coerce')

    # Filter for assemblies with contigN50 >= 100000 (Annotation Candidates)
    annotation_candidates = df_filtered_levels[df_filtered_levels['Contig N50'] >= 100000].groupby(
        'Year').size().reset_index(
        name='Annotation Candidates')

    # Extract annotation year from 'last_geneset_update'
    if 'last_geneset_update' in filtered_df.columns:
        filtered_df['Year'] = pd.to_datetime(filtered_df['last_geneset_update']).dt.year
        annotated_per_year = filtered_df.groupby('Year').size().reset_index(name='Number of Annotated Genomes')
    else:
        print("No 'last_geneset_update' column found in annotations.")
        return pd.DataFrame()

    # Extract Ensembl release year
    if 'ensembl_release_date' in filtered_df.columns:
        filtered_df['Ensembl Release Year'] = pd.to_datetime(filtered_df['ensembl_release_date']).dt.year
        ensembl_release_summary = filtered_df.groupby('Ensembl Release Year').size().reset_index(
            name='Number of Ensembl Releases')
    else:
        ensembl_release_summary = pd.DataFrame()

    # Merge all summaries
    yearly_summary = pd.merge(assemblies_per_year, annotated_per_year, on='Year', how='outer').fillna(0)
    yearly_summary = pd.merge(yearly_summary, annotation_candidates, on='Year', how='outer').fillna(0)
    yearly_summary = pd.merge(yearly_summary, ensembl_release_summary, left_on='Year', right_on='Ensembl Release Year', how='outer').fillna(0).drop(
        columns=['Ensembl Release Year'])

    return yearly_summary


def bin_by_assembly_level(df):
    """Bins assemblies based on their assembly level."""
    if 'Assembly level' not in df.columns:
        print("Column 'Assembly level' not found in dataset.")
        return pd.DataFrame()  # Return empty DataFrame if column is missing

    # Group by Assembly level and count occurrences
    level_summary = df.groupby('Assembly level', observed=False).size().reset_index(name='Number of Assemblies')

    return level_summary



def main():
    parser = argparse.ArgumentParser(description='Fetches annotations and assemblies')
    parser.add_argument("--release_date", type=str, default="2019-01-01", help="Release date in YYYY-MM-DD format (default: 2019-01-01).")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save output files.")
    parser.add_argument("--taxon_id", type=str, required=False, help="Taxon ID for filtering (e.g., 40674 for Mammalia).")
    parser.add_argument('--asm_level', type=str, nargs='+', help="Assembly level options: 'Contig', 'Scaffold', 'Chromosome', 'Complete genome'.")
    parser.add_argument('--asm_type', type=str, nargs='+', help="Assembly type: 'haploid', 'alternate-pseudohaplotype', 'unresolved-diploid', 'haploid-with-alt-loci', 'diploid'.")
    parser.add_argument('--contig_n50', type=float, help="Contig N50")

    args = parser.parse_args()

    release_date = args.release_date
    output_dir = args.output_dir
    taxon_id = args.taxon_id
    asm_level = args.asm_level
    asm_type = args.asm_type

    os.makedirs(output_dir, exist_ok=True)

    # Define metrics to consider
    metric_thresholds = {k: v[0] if isinstance(v, list) else v for k, v in vars(args).items() if
                         v is not None and k not in ["bioproject_id", "output_dir", "asm_level", "asm_type",
                                                     "release_date"]}

    all_metrics = ["gc_percent", "total_sequence_length", "contig_n50", "number_of_contigs", "number_of_scaffolds",
                   "scaffold_n50", "genome_coverage"]

    # Fetch assemblies released in the past 5 years
    df_wide, summary_df, df_info_result, df_gca_list = get_filtered_assemblies(None, metric_thresholds, all_metrics,asm_level, asm_type, release_date, taxon_id)


    if isinstance(df_wide, str):  # Check if there's an error message
        print(df_wide)
        return

    # Bin assemblies by Assembly level
    level_summary = bin_by_assembly_level(df_wide)

    # Fetch dataset annotations using the same release date
    filtered_df=get_annotation(release_date)

    #Add FTP links to dataframe
    #filtered_df['FTP'] = filtered_df['GCA'].apply(
        #lambda gca: f"https://ftp.ebi.ac.uk/pub/ensemblorganisms/{filtered_df.loc[filtered_df['GCA'] == gca, 'scientific_name'].values[0]}/{gca}/")

    # Bin the filtered results by genebuild.method
    method_summary = bin_by_genebuild_method(filtered_df)

    yearly_summary = generate_yearly_summary(df_wide, df_info_result, filtered_df)

    #Check annotation status
    # Add an 'Annotated' column to df_info_result
    df_info_result['Annotated'] = df_info_result['GCA'].isin(filtered_df['GCA']).map({True: 'Yes', False: 'No'})

    print("\nDataset filtered:")
    print(method_summary)

    print("\nDataset yearly summary:")
    print(yearly_summary)

    print("\nAnnotated:")
    print(filtered_df)

    print("\nAssembly Level Summary:")
    print(level_summary)


# Save results to CSV files
    method_summary.to_csv(os.path.join(output_dir, "annotation.csv"), index=False)
    yearly_summary.to_csv(os.path.join(output_dir, "yearly_summary.csv"), index=False)
    filtered_df.to_csv(os.path.join(output_dir, "filtered_df.csv"), index=False)
    df_wide.to_csv(os.path.join(output_dir, "assemblies.csv"))
    df_info_result.to_csv(os.path.join(output_dir, "assembly_annotation_status.csv"))
    level_summary.to_csv(os.path.join(output_dir, "assembly_level_summary.csv"), index=False)

if __name__ == "__main__":
    main()
