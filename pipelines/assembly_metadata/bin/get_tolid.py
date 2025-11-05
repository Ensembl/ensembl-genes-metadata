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

import requests
import argparse
import pymysql
import json
from tenacity import retry, stop_after_attempt, wait_random
import logging
import os

def connect_db(query:str, metadata_params) -> tuple:
    """
    Connect to the assembly metadata database and execute a query

    Args:
        query (str): query to be executed

    Returns:
        tuple: results from the query's execution
    """
    logging.info(f"Querying metadata database with: {query}")
    con = pymysql.connect(**metadata_params)
    with con.cursor() as cursor:
        cursor.execute(query)
        data = cursor.fetchone()
    con.close()
    return data

@retry(stop=stop_after_attempt(5), wait=wait_random(min=1, max=10))
def connection_TolID(uri: str) -> requests.Response:
    response = requests.get(uri)
    response.raise_for_status()
    return response

def get_tolid(taxon: int) -> str:
    """Retrieves the ToLID prefix for a given taxon from the ToLID API

    Args:
        taxon (int): lowest taxon id of the assembly

    Returns:
        str: ToLID prefix for the given taxon
    """
    
    logging.info(f"Connecting to ToLID API to retrieve prefix for taxon {taxon}")
    uri = f"https://id.tol.sanger.ac.uk/api/v2/species?taxonomyId={taxon}"
    response = connection_TolID(uri)
    data = response.json()
    logging.info(f"Response from ToLID API: {data}")
    
    if data['totalNumSpecies']==0:
        logging.info(f"No species found for taxon {taxon}")
        tolid_prefix = ""
    else:
        tolid_prefix = data['species'][0]['prefix']
    
    return tolid_prefix

def main():
    """Module's entry point
    """
    
    logging.basicConfig(filename="get_tolid.log", level=logging.DEBUG, filemode='w',
                        format="%(asctime)s:%(levelname)s:%(message)s")
    
    parser = argparse.ArgumentParser(prog='get_tolid.py',
        description='Get ToLID prefix for a given assembly based on the taxon ID')
    parser.add_argument('--accession', 
                        type=str,
                        help='Full GCA assembly accession')
    
    parser.add_argument('--metadata',
                        help='JSON file with metadata database connection parameters')
    
    args = parser.parse_args()
    logging.info(args)
    
    if args.metadata:
        if os.path.exists(args.metadata):
            with open(args.metadata, 'r') as file:
                metadata_params = json.load(file)
        else:
            raise ValueError("Metadata params json file does not exist")
    
    logging.info(f"Connecting to metadata database to retrieve lowest_taxon_id and assembly_id for {args.accession}")
    gca_chain = args.accession.split('.')[0]
    gca_version = args.accession.split('.')[1]
    query_taxon =f'SELECT assembly.assembly_id, species.species_taxon_id FROM assembly INNER JOIN species ON species.lowest_taxon_id = assembly.lowest_taxon_id WHERE gca_chain = "{gca_chain}" AND gca_version = "{gca_version}"'
    assembly_id, taxon = connect_db(query_taxon, metadata_params)
    
    tolid = get_tolid(taxon)
    
    logging.info(f"Connecting to metadata database to retrieve organism_id for {args.accession}")
    query_organism = f'SELECT organism_id FROM organism WHERE assembly_id = "{assembly_id}"'
    organism_id = connect_db(query_organism, metadata_params)
    
    dtol_json = {
        'organism' : {
            'organism_id': organism_id[0], 
            'tol_prefix': tolid
            }
        }
    
    logging.info(f"Saving tolid json for {args.accession}")
    with open(f'{args.accession}_tolid.json', 'w') as file:
        json.dump(dtol_json, file)
    file.close()
    
if __name__ == '__main__':
    main()