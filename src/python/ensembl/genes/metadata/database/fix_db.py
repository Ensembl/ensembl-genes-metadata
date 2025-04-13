import pymysql
import pandas as pd
from tqdm import tqdm
import argparse
import re
from typing import Dict
import json
import requests
from tenacity import (
    retry,
    stop_after_attempt,
    wait_exponential,
    wait_random,
)
BATCH_SIZE = 50 

def connect_to_db(host, user, password, database, port=3306):
    """
    Establishes a connection to the MySQL database.

    :param host: The hostname or IP address of the MySQL server.
    :param user: The username to connect to the database.
    :param password: The password for the user.
    :param db: The name of the database to connect to.
    :param port: The port of the MySQL server (default is 3306).
    :return: A pymysql connection object.
    """
    try:
        connection = pymysql.connect(
            host=host,
            user=user,
            password=password,
            database=database,
            port=port
        )
        return connection
    except pymysql.MySQLError as e:
        print(f"Error connecting to MySQL: {e}")
        return None
def get_sample_info(accession: str) -> str:
    """Get info about sample name and description for the run accession"""
    biosample_url = f"https://www.ebi.ac.uk/biosamples/samples/{accession}"
    print(f"Fetching data from {biosample_url}")
    try:
        response = requests.get(biosample_url, timeout=10)
        response.raise_for_status()  # Raise an HTTPError if the request was not successful
        biosample_data = response.json()
        characteristics = biosample_data.get("characteristics", {})
        sample = ""
        if "tissue" in characteristics:
            sample = characteristics["tissue"][0]["text"].lower()
        elif "source_name" in characteristics:
            sample = characteristics["source_name"][0]["text"].lower()
        sample = re.sub(r"[ ;\(\)\/\\]", " ", sample)
        # remove punctuation
        sample = re.sub(r"[!\"#$%&()*\+,\-\'.\/:;<=>?@\[\]^`{|}~]", "", sample)
        return sample.strip()
    except (requests.RequestException, requests.HTTPError, ConnectionError, requests.Timeout) as e:
        print(f"An error occurred while fetching data from {biosample_url}: {str(e)}")
        # Handle the error here, you can log it or take other appropriate actions.
        return ''
def request_data(run_accession: str) -> str:
    """Make an HTTP request for the metadata of the given run_accession.

    Args:
        run_accession (str): Unique run accession.
        fields (list): List of fields of interest.

    Returns:
        str: TSV file containing metadata.
    """

    query = 'run_accession="' + run_accession + '"'

    data = {
        "result": "read_run",
        "query": query,
        "fields": "library_layout",
        "format": "tsv",
    }
    # @retry(stop=stop_after_attempt(3), wait=wait_fixed(5))
    @retry(
        stop=stop_after_attempt(5),
        wait=wait_exponential(multiplier=1, min=5, max=60) + wait_random(min=0, max=5),
    )
    def make_request():
        try:
            response = requests.post("https://www.ebi.ac.uk/ena/portal/api/search", data=data, timeout=20)
            response.raise_for_status()
            print(f"Request successful for {run_accession}")
            table = response.text.split("\n")
            meta_data = table[1].split("\t")
            output_data = {}
            fields = ["run_accession","library_layout"]
            for key, value in zip(fields, meta_data):
                output_data[key] = value
            print(f"Meta data: {meta_data}")
            paired= 1 if output_data["library_layout"] == "PAIRED" else 0
            print(f"Paired: {paired}")
            return paired

        except requests.RequestException as e:
            raise RuntimeError(f"Error fetching metadata: {e}") from e

    try:
        result = make_request()
        return result
        # Further processing of result can be done here
    except RuntimeError as e:
        print("Error retriving data")
        return str(e)

if __name__ == "__main__":
    # Connection details
    db_config = {
        'host': 'mysql-ens-genebuild-prod-1',
        'user': 'ensadmin',
        'password': 'ensembl',
        'database': 'gb_transcriptomic_registry_copy_110325',
        "port": 4527
    }
    # Load the rows to fix (e.g., 150k rows)
    query = "SELECT run_accession,sample_accession  FROM run"
    # Connect to the database
    db_connection = connect_to_db(**db_config)
        # Clean the input text
    cursor = db_connection.cursor()
    if db_connection:
        df = pd.read_sql(query, db_connection)
    print(f"Loaded {len(df)} rows.")
    # Update back to DB
    update_query_paired = "UPDATE run set paired = %s where run_accession = %s"
    update_query_biosample = "UPDATE run set sample_tissue = %s where sample_accession = %s"
    counter = 0
    for _, row in tqdm(df.iterrows(), total=len(df)):
        try:
            run_acc = row["run_accession"]
            sample_acc = row["sample_accession"]
            print(f"Updating row {row['run_accession']}")
            paired_val = request_data(run_acc)
            tissue_val = get_sample_info(sample_acc)
            cursor.execute(update_query_paired, (paired_val, run_acc))
            cursor.execute(update_query_biosample, (tissue_val, sample_acc))
            counter += 1
            # Commit every BATCH_SIZE rows
            if counter % BATCH_SIZE == 0:
                db_connection.commit()
                print(f"Committed {counter} rows...")
            
            #cursor.execute(update_query_paired, (request_data(row["run_accession"]), row["run_accession"]))
            #cursor.execute(update_query_biosample, (get_sample_info(row["sample_accession"]), row["sample_accession"]))
            
        except Exception as e:
            print(f"Failed to update row {row['run_accession']}: {e}")
            
    # Commit changes
    db_connection.commit()
    cursor.close()
    db_connection.close()
    print("✅ Update complete.")
