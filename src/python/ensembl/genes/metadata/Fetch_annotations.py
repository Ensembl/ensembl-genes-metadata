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
    # Convert release_date to datetime and extract year in df_info_result
    df_info_result['Year'] = pd.to_datetime(df_info_result['Release date']).dt.year

    # Merge Year information into df_wide based on assembly_id
    df = df.merge(df_info_result[['GCA', 'Year']], on='GCA', how='left')

    # Count assemblies above contig level (scaffold, chromosome, complete genome)
    valid_levels = ["Scaffold", "Chromosome", "Complete genome"]
    assemblies_per_year = df[df['Assembly level'].isin(valid_levels)].groupby('Year').size().reset_index(
        name='Number of Assemblies')

    # Convert annotation release_date to datetime and extract year if column exists
    if 'last_geneset_update' in filtered_df.columns:
        filtered_df['Year'] = pd.to_datetime(filtered_df['last_geneset_update']).dt.year
        annotated_per_year = filtered_df.groupby('Year').size().reset_index(name='Number of Annotated Genomes')
    else:
        print("No 'last_geneset_update' column found in annotations.")
        return pd.DataFrame()

    # Merge both summaries
    yearly_summary = pd.merge(assemblies_per_year, annotated_per_year, on='Year', how='outer').fillna(0)

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
    parser.add_argument("--release_date", type=str, required=True, help="Release date in YYYY-MM-DD format.")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save output files.")

    args = parser.parse_args()

    release_date = args.release_date
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    # Define metrics to consider
    metric_thresholds = {}  # Empty, so no thresholds are applied
    all_metrics = ["gc_percent", "total_sequence_length", "contig_n50", "number_of_contigs", "number_of_scaffolds",
                   "scaffold_n50", "genome_coverage"]

    # Fetch assemblies released in the past 5 years
    df, summary_df, df_info_result, df_gca_list = get_filtered_assemblies(None, metric_thresholds, all_metrics, None, None, release_date)


    if isinstance(df, str):  # Check if there's an error message
        print(df)
        return

    # Bin assemblies by Assembly level
    level_summary = bin_by_assembly_level(df)


    print("\nAssembly Level Summary:")
    print(level_summary)

    # Save to CSV
    level_summary.to_csv(os.path.join(output_dir, "assembly_level_summary.csv"), index=False)


    # Fetch dataset annotations using the same release date

    filtered_df=get_annotation(release_date)

    # Bin the filtered results by genebuild.method
    method_summary = bin_by_genebuild_method(filtered_df)

    yearly_summary = generate_yearly_summary(df, df_info_result, filtered_df)


    print("\nDataset filtered:")
    print(method_summary)

    print("\nDataset yearly summary:")
    print(yearly_summary)

    print("\nAnnotated:")
    print(filtered_df)

# Save results to CSV files
    method_summary.to_csv(os.path.join(output_dir, "annotation.csv"), index=False)
    yearly_summary.to_csv(os.path.join(output_dir, "yearly_summary.csv"), index=False)
    filtered_df.to_csv(os.path.join(output_dir, "filtered_df.csv"), index=False)
    df.to_csv(os.path.join(output_dir, "df.csv"))
    df_info_result.to_csv(os.path.join(output_dir, "df_info_result.csv"))


if __name__ == "__main__":
    main()
