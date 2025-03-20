#  See the NOTICE file distributed with this work for additional information
#  regarding copyright ownership.
#
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#  http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
import logging
import argparse
import pymysql #type: ignore
import requests #type: ignore
from datetime import datetime
import json
from tenacity import retry, stop_after_attempt, wait_random

def execute_query(query, db_params):
    conn = pymysql.connect(**db_params)
    cursor = conn.cursor()
    cursor.execute(query)
    result = cursor.fetchall()
    cursor.close()
    conn.close()
    return result

@retry(stop=stop_after_attempt(5), wait=wait_random(min=1, max=10))
def connection_ncbi(uri: str) -> requests.Response:
    """Connects to NCBI API and retrieves data for a given URI
    Args:
        uri (str): URI to request data from NCBI API
    Returns:
        requests.Response: response object from NCBI API
    """
    response = requests.get(uri, timeout=60)
    response.raise_for_status()
    return response

def fetch_gca_list(metadata_params, full_screen):
    """ Get a list of GCAs to check their status and generate the update queries
    Args:
        metadata_params (dict): dictionary containing the connecting params for the gb_assembly_metadata database
        full_screen (boolean): if true it will update all the available records in the DB. Otherwise, it will update vertebrate or assemblies from the biodiversity projects 
    Returns:
        list: GCAs list
    """
    
    if full_screen:
        logging.info("Running in full screen mode")
        query_get_gca_list = """
        SELECT CONCAT(a.GCA_chain, '.', a.GCA_version) AS GCA
        FROM assembly a
        WHERE a.is_current = 'current';
        """
    else:
        logging.info("Running only high priority assemblies (vertebrates or relevant bioprojects)")
        query_get_gca_list = """
        SELECT DISTINCT CONCAT(a.GCA_chain, '.', a.GCA_version) AS GCA
        FROM bioproject b
        JOIN assembly a ON a.assembly_id = b.assembly_id
        JOIN species s ON a.lowest_taxon_id = s.lowest_taxon_id
        JOIN taxonomy t ON t.lowest_taxon_id = a.lowest_taxon_id
        WHERE s.scientific_name NOT LIKE '% sp.%'
            AND a.is_current = 'current'
            AND (t.taxon_class_id = '7711'
                OR b.bioproject_id IN (
                    'PRJNA533106', 'PRJEB40665', 'PRJEB61747',
                    'PRJEB43510', 'PRJEB47820', 'PRJNA813333', 'PRJNA489243'
        ));
        """
        
    #Get list from assembly_metadata and parse
    fetch_gca = execute_query(query_get_gca_list, metadata_params)
    gca_list = [gca[0] for gca in fetch_gca]
    logging.info(f"Total number of assemblies to check: {len(gca_list)}")
    
    return gca_list
    
def update_status(gca_list, metadata_params, ncbi_url):
    """Update the status of a given GCA list
    Args:
        gca_list (list): list of GCA to update
        metadata_params (dict): dictionary containing connecting parameters for gb_assembly_metadata
        ncbi_url (str): current NCBI API URL
    Returns:
        ist[list, list, list]: returns three list of GCAs that were modified: not_current_list, warnings_list, deleted_list
    """
    # List to store GCAs that will be updated
    not_current_list = []
    warnings_list = []
    deleted_list = []
    
    for accession in gca_list:
        uri = f"{ncbi_url}/genome/accession/{accession}/dataset_report?filters.exclude_atypical=false&filters.assembly_version=all_assemblies"
        #logging.info("URI: %s",uri)
        response = connection_ncbi(uri)
        data = response.json()
        if data == {}:
            logging.info(f"Assembly {accession} is deleted from NCBI")
            query_update_status = f"UPDATE assembly a SET is_current = 'suppressed' WHERE CONCAT(a.gca_chain, '.', a.gca_version) = '{accession}';"
            execute_query(query_update_status, metadata_params)
            deleted_list.append(accession)
        else:
            is_current = data['reports'][0]['assembly_info']['assembly_status']
            warning = data['reports'][0].get('assembly_info').get('atypical', {}).get('warnings', "NA")

            if is_current != "current":
                logging.info(f"Assembly {accession} is not current")
                query_update_status = f"UPDATE assembly a SET is_current = '{is_current}' WHERE CONCAT(a.gca_chain, '.', a.gca_version) = '{accession}';"
                execute_query(query_update_status, metadata_params)
                not_current_list.append(accession)
            if warning != "NA":
                logging.info(f"Assembly {accession} has warnings: {warning[0]}")
                query_update_status = f"UPDATE assembly a SET is_current = '{warning[0]}' WHERE CONCAT(a.gca_chain, '.', a.gca_version) = '{accession}';"
                execute_query(query_update_status, metadata_params)
                warnings_list.append(accession)
                
    return not_current_list, warnings_list, deleted_list

def main():
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    logging.basicConfig(
        filename=f"update_status_gca_records_{timestamp}.log",
        level=logging.DEBUG,
        filemode='w',
        format="%(asctime)s:%(levelname)s:%(message)s")

    parser = argparse.ArgumentParser(
        prog='update_status_gca_records.py',
        description='Update status of GCA records'
    )
    parser.add_argument('--metadata',
                        type=str,
                        help="Path to the metadata database params in json format")
    parser.add_argument('--ncbi_url',
                        type=str,
                        default="https://api.ncbi.nlm.nih.gov/datasets/v2",
                        help="NCBI API URL")
    parser.add_argument('--full_screen',
                        action='store_true',
                        help='It will check the assembly status of all available GCAs')
    args = parser.parse_args()
    logging.info(f"Script arguments: {args}")

    # Loading metadata DB params
    with open(args.metadata, 'r') as file:
        metadata_params = json.load(file)

    # Fetch list of genomes
    gca_list = fetch_gca_list(metadata_params, args.full_screen)

    # Update assembly status
    not_current_list, warnings_list, deleted_list = update_status(gca_list, metadata_params, args.ncbi_url)

    logging.info("Summary:")
    logging.info("Mode: only high priority assemblies (vertebrates or relevant bioprojects)")
    logging.info(f"Total number of assemblies to check: {len(gca_list)}")
    logging.info(f"Assemblies that were updated to previous or suppressed status: {len(not_current_list)}")
    logging.info(f"Assemblies with warnings: {len(warnings_list)}")
    logging.info(f"Assemblies that were deleted (not longer found on NCBI ): {len(deleted_list)}")

    # Saving list
    with open(f"list_not_current_{timestamp}.txt", "w") as file:
        for gca in not_current_list:
            file.write(gca + "\n")
    with open(f"list_warning_{timestamp}.txt", "w") as file:
        for gca in warnings_list:
            file.write(gca + "\n")
    with open(f"list_deleted_{timestamp}.txt", "w") as file:
        for gca in deleted_list:
            file.write(gca + "\n")     

if __name__ == "__main__":
    main()