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


"""The module will connect to NCBI API and Genebuild databases to determinate the list of 
GCA accession to register in the new assembly metadata database

Args:
    taxon (int): Taxon ID, by default it is 2759 (Eukaryote domain)
    file_path (str): path with the last database update (optional)

Returns:
    stdout: prints to the standard output a list of GCA accession
"""

import requests
import os
from datetime import datetime
import pymysql
import argparse
import logging
from tenacity import retry, stop_after_attempt, wait_random
from typing import Dict, List, Any, Tuple

from gb_db_params import METADATA_RPARAMS, REGISTRY_RPARAMS
from ncbi_params import NCBI_PARAMS

def set_date(taxon:int, ncbi_params: Dict[str, Any], file_path:str) -> Tuple[Dict[str, str], str]:
    """
    Set a date to retrieve new assemblies from NCBI. If a file with the last update date is provided, 
    date will be used to retrieve new assemblies. Otherwise, a default date (01/01/2019) will be used.

    Args:
        taxon (int): Taxon ID, by default it is 2759 (Eukaryote domain)
        ncbi_params (dict): NCBI API's parameters
        file_path (str): path with the last database update (optional) 

    Returns:
        dict: Updated NCBI API's parameters
        str: date to be used to retrieve assemblies 
    """
    logging.info('Checking date to retrieve assemblies')
    
    release_date = "01/01/2019"

    if file_path:
        with open(file_path, 'r') as file:
            for line in file:
                if int(line.split()[0]) == taxon:
                    release_date = line.split()[1]
                    ncbi_params['filters.first_release_date'] = release_date
                    logging.info(f"Last update record found! {line.split()[1]} will be used to retrieve new available assemblies for taxon {taxon}")
            else:
                logging.info(f"There is not last update record saved for taxon {taxon}. Using {release_date} as default date to retrieve assemblies")
    else:
        logging.info(f"Last update file not found. Using {release_date} as default date to retrieve assemblies")

    return ncbi_params, release_date

@retry(stop=stop_after_attempt(5), wait=wait_random(min=1, max=10))
def connection_ncbi(uri: str, params: Dict[str, str]) -> requests.Response:
    response = requests.get(uri, params=params)
    response.raise_for_status()
    return response

def fetch_gca_list(taxon: int, ncbi_params: Dict[str, str]) -> List[str]:
    """
    Fetch a list of GCA accession from NCBI API based on the taxon ID
    
    Args:
        taxon (int): Taxon ID, by default it is 2759 (Eukaryote domain)
        ncbi_params (dict): NCBI API's parameters

    Returns:
        list: list of GCA accessions
    """
    gca_list = [] 
    page_token=None
    uri = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon/{str(taxon)}/dataset_report"
    next_page = True
    
    while next_page:
        # Fetching data from NCBI API
        response = connection_ncbi(uri, ncbi_params)
        logging.info(f"URL: {response.url}")
        data = response.json()
        
        # Extracting GCA accessions
        gca_list.extend(report['accession'] for report in data['reports'])

        # Setting up the next page
        page_token = data.get('next_page_token')
        if page_token:
            ncbi_params['page_token'] = page_token
        elif page_token is None: # No more pages, break loop
            next_page = False
        else:
            raise ValueError("Page token is not being set correctly")  

    return set(gca_list)       

def build_db_query(release_date: str) -> Dict[str, str]:
    """ Build mysql the query to retrieve data from the db
    
    Args:
        release_date (str): date to be used to retrieve assemblies
    
    Returns: 
        dict: a dictionary containing a query for each database
    
    """
    
    release_date_sql = datetime.strptime(release_date, '%m/%d/%Y').strftime('%Y-%m-%d')
    
    query_registry =  """SELECT concat(assembly.chain, '.', assembly.version) AS GCA
    FROM assembly
    INNER JOIN meta
    ON assembly.assembly_id = meta.assembly_id
    WHERE meta.assembly_date >= DATE('{}') ; """.format(release_date_sql)
    
    query_metadata = """SELECT concat(gca_chain, '.', gca_version) 
    FROM assembly 
    WHERE release_date >= DATE('{}') ; """.format(release_date_sql)
    
    db_query = {
    'asm_metadata':
        {'db': 'gb_assembly_metadata_testing',
        'query': query_metadata},
    'asm_registry':
        {'db': 'gb_assembly_registry',
        'query': query_registry},
    }

    return db_query

def fetch_records_db(db_params: Dict[str, Any], query: str) -> List[str]:
    """Fetch assemblies that have been register after the last update

    Args:
        db_params (dict): database connection parameters
        query (str): mysql query to retrieve data 

    Returns:
        list: list of GCA recorded after the last update date
    """
    ## Connection to the DB
    with pymysql.connect(**db_params) as conn:
        with conn.cursor() as cur:
            # Querying database
            cur.execute(query)
            output= cur.fetchall() 
            # Processing results
            reg_gca = [row[0] for row in output]
            
    return reg_gca

def get_gca_register(db: str, db_query: Dict[str, Any], gca_list: List[str]) -> List[str]:
    """ Get list of GCA assemblies to register

    Args:
        db (str): option indicating what database use to retrieve accession and later compare
                    (asm_metadata, asm_registry, both)
        db_query (dict): database and its query
        gca_list (list): list of GCA accessions retrieved from NCBI API
    """
    accessions_records = []
    records_metadata = []
    records_registry = []
    
    if db in ['asm_metadata', 'both']:
        logging.info('Getting assemblies from assembly metadata database')
        records_metadata = fetch_records_db(METADATA_RPARAMS, query=db_query['asm_metadata']['query'] )
    elif db in ['asm_registry', 'both']:
        logging.info('Getting assemblies from assembly registry database')
        records_registry = fetch_records_db(REGISTRY_RPARAMS, query=db_query['asm_registry']['query'] )
        
    accessions_records = records_metadata + records_registry
    accessions_to_register = gca_list - set(accessions_records)
    
    return accessions_to_register

def main():
    """module's entry-point
    """
    
    logging.basicConfig(filename="fetch_new_assemblies.log", level=logging.DEBUG, filemode='w',
                    format="%(asctime)s:%(levelname)s:%(message)s")
            
    parser = argparse.ArgumentParser(prog='fetch_new_assemblies.py', 
                                    description='Identify new assemblies to register in the assembly metadata database. Assemblies are retrieved from NCBI API and compared with the assembly metadata and assembly registry databases.')
    
    parser.add_argument('--taxon', 
                        default=2759,
                        type=int,
                        help='Valid Taxon id: Eukaryota - 2759')
    parser.add_argument('--db',
                        default='both',
                        type=str,
                        help='Database to compare the assemblies (asm_metadata, asm_registry, both). Default is "both"')
    parser.add_argument('--file-path',
                        type=str,
                        help="File with the last update date")

    args = parser.parse_args()
    logging.info(args)
    
    if isinstance(args.taxon, int):
        taxon = args.taxon
    else:
        raise ValueError("Please enter a valid taxon numeric value.")
    
    if args.db not in ['asm_metadata', 'asm_registry', 'both']:
        raise ValueError("Please enter a valid database name (asm_metadata, asm_registry, both)")
    
    if args.file_path:
        if not os.path.exists(args.file_path):
            raise ValueError("Please enter a valid file path")
    else:
        logging.info("Default date will be used to retrieve assemblies")
    
    ncbi_params = NCBI_PARAMS
    
    ncbi_params, release_date = set_date(taxon, ncbi_params, args.file_path)
    gca_list = fetch_gca_list(taxon, ncbi_params)
    db_query = build_db_query(release_date)
    accessions_to_register = get_gca_register(args.db, db_query, gca_list)
    
    with open("assemblies_to_register.txt", 'w') as file:
        for accession in accessions_to_register:
            print(accession)
            file.write(accession + '\n')
    file.close()
        
    logging.info(f'Accessions to register: {len(accessions_to_register)}. Please note that some assemblies might belong to unspecified species')

if __name__ == '__main__':
    main()