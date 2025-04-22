import os
import json
import pymysql
import logging
import pandas as pd


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



def main(taxonomy_dict):
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
    print(trans_taxon_ids)

    trans_placeholders = ','.join(['%s'] * len(trans_taxon_ids))
    where_clause = f"WHERE r.taxon_id IN ({trans_placeholders})"

    query = f"""
        SELECT m.taxon_id, r.qc_status, m.last_check, 
        a.uniquely_mapped_reads_percentage AS star_uniquely_mapped_reads, a. percentage_reads_mapped_to_multiple_loci AS star_mapped_to_multiple_loci,
        a.percentage_reads_unmapped_too_short AS star_unmapped_reads,
        d.basic_statistics AS fastqc_basic, d.per_base_sequence_quality AS fastqc_per_base, d.total_sequences AS fastqc_total_seq, d.sequence_length AS fastqc_sequence_length,
        d.gc_content AS fastqc_gc_content, d.file_name AS fastqc_file_name
        FROM run r
        JOIN meta m ON m.taxon_id = r.taxon_id
        JOIN align a ON a.run_id = r.run_id
        JOIN data_files as d ON d.run_id = r.run_id
        {where_clause};
    """

    cursor.execute(query, trans_taxon_ids)
    trans_results = cursor.fetchall()

    trans_df = pd.DataFrame(trans_results)

    # Split the DataFrame into two based on _1 and _2 in fastqc_file_name
    trans_df_1 = trans_df[trans_df['fastqc_file_name'].str.endswith('_1')].copy()
    trans_df_2 = trans_df[trans_df['fastqc_file_name'].str.endswith('_2')].copy()

    # Add a suffix to all fastqc columns (except file name) to indicate _1 or _2
    fastqc_cols = [col for col in trans_df.columns if col.startswith('fastqc_') and col != 'fastqc_file_name']
    trans_df_1 = trans_df_1.rename(columns={col: f"{col}_1" for col in fastqc_cols})
    trans_df_2 = trans_df_2.rename(columns={col: f"{col}_2" for col in fastqc_cols})

    # Drop the fastqc_file_name column since it's now implied in the suffix
    trans_df_1.drop(columns=['fastqc_file_name'], inplace=True)
    trans_df_2.drop(columns=['fastqc_file_name'], inplace=True)

    # Merge the two on shared keys (e.g., taxon_id, run_id or whatever is common)
    # Assuming 'r.run_id' is the common key
    merged_df = pd.merge(trans_df_1, trans_df_2, on=['taxon_id', 'qc_status', 'last_check', 'star_uniquely_mapped_reads',
                                         'star_mapped_to_multiple_loci', 'star_unmapped_reads'], how='outer')


    found_ids = trans_df['taxon_id'].nunique()
    logging.info(
        f"Transcriptomic assessment data retrieved for {found_ids} taxon IDs out of {len(trans_taxon_ids)} provided.")

    missing_count = len(trans_taxon_ids) - found_ids
    if missing_count > 0:
        logging.info(f"{missing_count} taxon IDs had no transcriptomic assessment data.")

    return merged_df


if __name__ == "__main__":
    main(taxonomy_dict)