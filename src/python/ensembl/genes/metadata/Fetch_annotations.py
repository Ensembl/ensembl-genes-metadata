import pandas as pd
import json
import pymysql
import os
from GB_metadata_time import get_filtered_assemblies


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


def get_annotation():
    db_config_key = "prod"
    conn = connect_db(db_config_key)
    cursor = conn.cursor()

    query = f"""
                SELECT d.name, dss.name AS db_name
                FROM dataset d
                JOIN dataset_source dss ON dss.dataset_source_id = d.dataset_source_id
                WHERE d.name = 'genebuild'
            """

    cursor.execute(query)
    results = cursor.fetchall()
    conn.close()

    df_prod = pd.DataFrame(results)
    return df_prod


def get_core_db_info(release_date):
    """Retrieve information from metadata table in different databases."""
    # Get the db_name from get_annotation
    annotation_data = get_annotation()

    all_results = []  # List to store results from all databases

    # Loop through all rows in annotation_data (which contains different db_names)
    for index, row in annotation_data.iterrows():
        db_name = row['db_name']  # Get db_name from the DataFrame
        print(f"Checking database: {db_name}")

        # Only query for specific meta_keys: 'genebuild.last_geneset_update' and 'assembly.accession'
        query_st = f"""
                    SELECT meta_key, meta_value
                    FROM meta
                    WHERE meta_key IN ('genebuild.last_geneset_update', 'assembly.accession', 'genebuild.method')
                    """

        db_configs = ["st5", "st6"]  # Check both st5 and st6 configurations

        for config in db_configs:
            try:
                # Attempt to connect to the database
                conn = connect_db(config)
                with conn.cursor() as cur:
                    # Check if the database exists before switching to it
                    cur.execute(f"SHOW DATABASES LIKE '{db_name}'")
                    result = cur.fetchone()

                    if result:
                        cur.execute(f"USE {db_name}")  # Switch to the relevant database
                        cur.execute(query_st)
                        output = cur.fetchall()

                        if output:
                            # Convert the results into a pandas DataFrame
                            df_output = pd.DataFrame(output, columns=['meta_key', 'meta_value'])

                            # Pivot the DataFrame to have meta_keys as columns
                            df_pivoted = df_output.pivot_table(columns='meta_key', values='meta_value', aggfunc='last').reset_index()
                            # Reset index for cleaner output
                            df_pivoted.reset_index(drop=True, inplace=True)
                            print(df_pivoted.head())

                            # Check if 'genebuild.last_geneset_update' exists before attempting operations on it
                            if 'genebuild.last_geneset_update' in df_pivoted.columns:
                                # Ensure genebuild.last_geneset_update has the correct format by adding '-01' for day precision
                                df_pivoted['genebuild.last_geneset_update'] = pd.to_datetime(
                                    df_pivoted['genebuild.last_geneset_update'] + '-01', errors='coerce')
                                # Debug: Print before filtering
                                print("\nBefore filtering:")
                                print(df_pivoted[['genebuild.last_geneset_update']].dropna().head())

                                # Convert the release_date to a datetime object
                                release_date = pd.to_datetime(release_date)

                                # Apply the release_date filter to 'genebuild.last_geneset_update'
                                df_filtered = df_pivoted[df_pivoted['genebuild.last_geneset_update'] >= release_date]

                                # Debug: Print after filtering
                                print("\nAfter filtering:")
                                print(df_filtered[['genebuild.last_geneset_update']].head())

                                # Log the number of retained rows
                                print(f"\nNumber of assemblies after filtering: {df_filtered.shape[0]}")

                                # Append this database's filtered data to all_results
                                all_results.append(df_filtered)
                                print(f"Data fetched for database: {db_name}")
                            else:
                                print(f"Column 'genebuild.last_geneset_update' not found in database: {db_name}")

                        else:
                            print(f"No data found in database: {db_name}")
                    else:
                        print(f"Database {db_name} does not exist in {config} configuration.")
            except pymysql.MySQLError as e:
                print(f"Database connection failed for config {config}: {e}")
            finally:
                if 'conn' in locals() and conn.open:
                    conn.close()

    # Check if we found any results for all databases
    if all_results:
        # Concatenate all the individual results into one DataFrame
        df_all_filtered = pd.concat(all_results, ignore_index=True)
        return df_all_filtered
    else:
        print("No data found in any of the databases.")
        return pd.DataFrame()  # Return an empty DataFrame if no results were found



def bin_by_contig_n50(df, bin_size=5000):
    """Bins assemblies based on contig N50 value."""
    bins = [0] + list(range(bin_size, df['Contig N50'].max() + bin_size, bin_size))
    labels = [f"{bins[i]}-{bins[i + 1] - 1}" for i in range(len(bins) - 1)]
    df['Contig N50 Bin'] = pd.cut(df['Contig N50'], bins=bins, labels=labels, right=False)
    return df


def summarize_assemblies(df):
    """Summarize the data by assembly level and contig N50 bins."""
    # Bin by assembly level
    level_summary = df.groupby('Assembly level',observed=False).size().reset_index(name='Number of Assemblies')

    # Bin by Contig N50
    contig_n50_summary = df.groupby('Contig N50 Bin',observed=False).size().reset_index(name='Number of Assemblies')

    # Other useful statistics (e.g., avg contig N50)
    avg_contig_n50 = df.groupby('Assembly level',observed=False)['Contig N50'].mean().reset_index(name='Avg Contig N50')

    return level_summary, contig_n50_summary, avg_contig_n50


def bin_by_genebuild_method(df):
    """Bins assemblies based on genebuild.method."""
    if 'genebuild.method' not in df.columns:
        print("Column 'genebuild.method' not found in dataset.")
        return pd.DataFrame()  # Return empty DataFrame if column is missing

    # Group by genebuild.method and count occurrences
    method_summary = df.groupby('genebuild.method', observed=False).size().reset_index(name='Number of Assemblies')

    return method_summary


def main():
    release_date = "2019-01-01"

    # Define metrics to consider
    metric_thresholds = {}  # Empty, so no thresholds are applied
    all_metrics = ["gc_percent", "total_sequence_length", "contig_n50", "number_of_contigs", "number_of_scaffolds",
                   "scaffold_n50", "genome_coverage"]

    # Fetch assemblies released in the past 5 years
    df, refseq_count, summary_df, info_result, df_gca_list = get_filtered_assemblies(release_date, metric_thresholds,
                                                                                     all_metrics, None, None)


    if isinstance(df, str):  # Check if there's an error message
        print(df)
        return

    # Bin assemblies by Contig N50
    df = bin_by_contig_n50(df)

    # Summarize the data
    level_summary, contig_n50_summary, avg_contig_n50 = summarize_assemblies(df)

    # Fetch dataset annotations using the same release date

    filtered_df=get_core_db_info(release_date)

    # Bin the filtered results by genebuild.method
    method_summary = bin_by_genebuild_method(filtered_df)

    # Print the results
    print("Summary by Assembly Level:")
    print(level_summary)

    print("\nSummary by Contig N50 Bin:")
    print(contig_n50_summary)

    print("\nAverage Contig N50 by Assembly Level:")
    print(avg_contig_n50)

    print("\nDataset filtered:")
    print(method_summary)

# Optionally, save results to CSV files
    level_summary.to_csv('/Users/lazar/Documents/Requests/Leanne/Assemblies_last_5_years/assembly_level_summary.csv', index=False)
    contig_n50_summary.to_csv('/Users/lazar/Documents/Requests/Leanne/Assemblies_last_5_years/contig_n50_summary.csv', index=False)
    avg_contig_n50.to_csv("/Users/lazar/Documents/Requests/Leanne/Assemblies_last_5_years/avg_contig_n50_by_level.csv", index=False)
    method_summary.to_csv("/Users/lazar/Documents/Requests/Leanne/Assemblies_last_5_years/annotation.csv", index=False)


if __name__ == "__main__":
    main()
