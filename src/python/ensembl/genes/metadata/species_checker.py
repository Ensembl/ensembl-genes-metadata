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

"""
The module get the additional information for the species table. It will add the species_taxon_id, parlance_name 
and species_prefix keys to the json-like (.tmp) species file.

Raises:
    ValueError: invalid taxonomy rank
    ValueError: multiple prefix detected

Returns:
    str: a json file with the species information
"""

import json
import argparse
import logging
import pymysql #type: ignore
import requests #type: ignore
import string
import random
import os
from tenacity import retry, stop_after_attempt, wait_random
from typing import Optional

@retry(stop=stop_after_attempt(10), wait=wait_random(min=1, max=20))
def connection_ncbi(uri: str) -> requests.Response:
    response = requests.get(uri)
    response.raise_for_status()
    return response

def get_taxon_data(taxon_id:str, ncbi_url) -> dict:
    """It connects to the NCBI API to retrieve the taxonomy data of the lowest taxon id.

    Args:
        taxon_id (int): lowest taxon id of the assembly
        ncbi_url (str): valid url to connect to NCBI API, taxonomy endpoint

    Returns:
        dict: response from NCBI API
    """
    uri = f"{ncbi_url}/taxonomy/taxon/{taxon_id}/dataset_report"
    logging.info(f"URI: {uri}")
    response = connection_ncbi(uri)
    taxon_data = response.json()

    return taxon_data

def get_taxon_classification(taxon_data) -> dict:
    """
    This function retrieves the taxonomy classification of the lowest taxon id.
    It returns the classification dictionary that will be user to update the taxonomy table.

    Args:
        lowest_taxon_id (str): the taxonomy id of the assembly, it is obtained from the assembly NCBI report

    """
    classification = taxon_data['reports'][0]['taxonomy']['classification']

    classification_dic = {}
    for rank in classification:
        if rank not in ['domain', 'superkingdom']:
            classification_dic.update({classification[rank]['id']: rank})

    return classification_dic

def species_taxon(taxon_data, taxon_id) -> tuple[int, bool]:
    """
    This functions checks if the lowest taxon id is a species or a infraspecific taxon rank.
    It returns the species taxon id to be used in the species table.

    Args:
        lowest_taxon_id (str): the taxonomy id of the assembly, it is obtained from the assembly NCBI report

    Raises:
        ValueError: when the provided value is not a valid taxonomy rank. Valid values are species taxon ID or infraspecific taxon ID

    Returns:
        str: species taxon id, it could be the same lowest taxon id value
    """

    taxon_exists = True
    try:
        taxonomy = taxon_data['reports'][0]['taxonomy']['rank']
        if taxonomy in ['SUBSPECIES', 'STRAIN', 'VARIETAS', 'GENOTYPE', 'ISOLATE', 'FORMA', 'FORMA_SPECIALIS']:
            species_taxon_id = taxon_data['reports'][0]['taxonomy']['classification']['species']['id']
            logging.info("The assembly is a infraspecific taxon %s", taxonomy)
        elif taxonomy == 'SPECIES':
            species_taxon_id = taxon_id
            logging.info("The assembly is a species taxon rank ")
        else:
            raise ValueError(f"Incorrect taxonomy ({taxonomy})")
    except KeyError:
        if 'errors' in taxon_data['reports'][0]:
            species_taxon_id = 0
            taxon_exists = False
            logging.info("Taxon %s is not a recognized NCBI Taxonomy name", taxon_id)
        elif 'taxonomy' in taxon_data['reports'][0]:
            logging.info('Taxon do not have Rank available, retrieving information from another section of the report')
            species_taxon_id = taxon_data['reports'][0]['taxonomy']['classification']['species']['id']
        else:
            raise KeyError("Taxon %s retrieves an unexpected report", taxonomy)

    return species_taxon_id, taxon_exists

def get_parlance_name(sci_name: str, enscode) -> str:
    """
    Search in the snp_static.txt file from the core_meta_update repository the parlance name 
    for the scientific name provided

    Args:
        sci_name (str): scientific name of the species
        enscode (str): ENSCODE variable path

    Returns:
        str: parlance name when available
    """

    parlance_file = f"{enscode}/ensembl-genes/src/python/ensembl/genes/metadata/snp_static.txt"
    data_dict = {}

    logging.info("Reading parlance name file (snp_static.txt) to look for a match")
    with open(parlance_file, 'r') as file:
        for line in file:
            key, value = line.rsplit("\t", 1)
            data_dict[key.strip()] = value.strip()

    parlance_name = data_dict.get(sci_name, "")

    return parlance_name

def create_prefix(registy_params, metadata_params) -> str:
    """
    It creates a unique species prefix. The prefix is the particle ENS plus a combination of 3 (or 4) random letters.
    To ensure the uniqueness of the prefix, it checks the existing prefixes in the registry and metadata databases.

    Returns:
        str: unique species prefix
    """
    
    # Getting existing prefix from registry db
    conn = pymysql.connect(**registy_params)
    cur  = conn.cursor()
    cur.execute("SELECT DISTINCT species_prefix FROM assembly")
    prefix_registry = cur.fetchall() 
    cur.close()
    
    # Getting existing prefix from metadata db
    conn = pymysql.connect(**metadata_params)
    cur  = conn.cursor()
    cur.execute("SELECT DISTINCT species_prefix FROM species")
    prefix_metadata = cur.fetchall() 
    cur.close()
    
    prefix_list = [item[0] for item in (prefix_registry + prefix_metadata) ]
    existing_prefix  = list(set(prefix_list))
    
    prefix = ''
    
    while not prefix or prefix in existing_prefix:
        letters = string.ascii_uppercase
        prefix = 'ENS'+ ''.join([random.choice(letters) for i in range(3)])
        
        if len(existing_prefix)>=26^3:
            prefix = 'ENS'+ ''.join([random.choice(letters) for i in range(4)])
    
    return prefix

def get_species_prefix(taxon_id:str, registy_params, metadata_params) -> Optional[str]:
    """
    This function retrieves the species prefix from the assembly registry and metadata databases.
    If the prefix is not found, it creates a new one. There are special cases where the prefix is predefined.
    - Canis lupus (wolf) -> ENSCAF
    - Canis lupus familiaris (Domestic dog) -> ENSCAF
    - Heterocephalus glaber (naked mole rat) -> ENSHGL

    Args:
        taxon_id (str): lowest taxon id

    Raises:
        ValueError: _description_

    Returns:
        str: unique species prefix that already exist or a new one when no prefix is found.
    """

    # Special cases
    special_cases = {
        '9612': 'ENSCAF', # Canis lupus (wolf)
        '9615' : 'ENSCAF', # Canis lupus familiaris (Domestic dog)
        '10181': 'ENSHGL' #  Heterocephalus glaber (naked mole rat)
        }

    if str(taxon_id) in special_cases:
        logging.info("The prefix is a special case")
        species_prefix = special_cases.get(str(taxon_id))
    else:

        logging.info(f"search prefix or create a new one for {taxon_id}")

        # Search prefix in DB (gb_assembly_registry)
        conn = pymysql.connect(**registy_params)
        cur  = conn.cursor()
        query = f"SELECT DISTINCT species_prefix FROM assembly WHERE taxonomy = {taxon_id}"
        cur.execute(query)
        output_registry = cur.fetchall()
        cur.close()

        # Search prefix in DB (gb_assembly_metadata)
        conn = pymysql.connect(**metadata_params)
        cur  = conn.cursor()
        query = f"SELECT DISTINCT species_prefix FROM species WHERE lowest_taxon_id = {taxon_id}"
        cur.execute(query)
        output_metadata = cur.fetchall()
        cur.close()

        # Combine output and get list of unique values
        output = [item[0] for item in (output_registry + output_metadata) ]
        prefix_list = list(set(output))

        # no prefix, create new prefix
        if len(prefix_list) == 0:
            logging.info(f"Getting a new prefix for taxon id: {taxon_id}")
            species_prefix = create_prefix(registy_params, metadata_params)
        # unique prefix detected
        elif len(prefix_list) == 1:
            logging.info(f"Unique prefix detected for taxon id: {taxon_id}")
            species_prefix = prefix_list[0]
        # Multiple prefix detected, check species.
        else:
            raise ValueError(f"The taxon {taxon_id} is already registered but multiple prefix were detected: {prefix_list}")

    return species_prefix

def main():
    """Module's entry point.
    """
    logging.basicConfig(filename="species_checker.log", level=logging.DEBUG,
                        filemode='w', format="%(asctime)s:%(levelname)s:%(message)s")
    parser = argparse.ArgumentParser(prog="species_checker.py",
                                    description=
                                    """
                                    Add species_taxon_id and parlance_name to the species main key to json-like (.tmp) species file.
                                    It all retrieve the taxonomy classification to be store in the taxonomy table.
                                    The module when used in the pipeline requires --json-path, --ncbi_url and --enscode.
                                    Additionally, to be used as standalone script, it requires --taxonomy_update and --taxon_id.
                                    """)
    parser.add_argument("--json-path",
                        type=str,
                        help="Path to the JSON-like (.tmp) species file")
    parser.add_argument('--ncbi_url',
                        type=str,
                        help='NCBI API URL')
    parser.add_argument('--enscode',
                        type=str,
                        help='ENSCODE path' )
    parser.add_argument('--taxon_id',
                        type=int,
                        help='Lowest taxon id of the species ')
    parser.add_argument('--taxonomy_update',
                        action='store_true',
                        help='Update taxonomy table')

    args = parser.parse_args()

    # Loading NCBI API URL
    if not args.ncbi_url:
        raise ValueError("Please enter a valid URL for NCBI API")

    logging.info(f"Loading file: {args.json_path}")

    # Loading Species dictionary
    if not args.taxonomy_update:
        with open(args.json_path, 'r') as file:
            species_dict = json.load(file)
        file.close()

    if not args.enscode and not args.taxonomy_update:
        raise ValueError("Please enter a valid path for ENSCODE")

    # Retrieve missing information to add keys to species json
    if not args.taxonomy_update:
        logging.info(f"Getting key values for the species: {species_dict['species']['scientific_name']}")
        # Get taxon data from NCBI API
        taxon_data = get_taxon_data(species_dict['species']['lowest_taxon_id'], args.ncbi_url)
        species_taxon_id, taxon_exists = species_taxon(taxon_data, species_dict['species']['lowest_taxon_id'])
        if taxon_exists:
            parlance_name = get_parlance_name(species_dict['species']['scientific_name'], args.enscode)
            species_prefix = ""
            taxon_classification = get_taxon_classification(taxon_data)
            taxon_classification_check=True
            #species_prefix = get_species_prefix(species_dict['species']['lowest_taxon_id'], registy_params, metadata_params)
        else:
            logging.info("Taxon do not exist in taxonomy: invalid lowest taxon id or assembly should be suppressed")
            logging.info("Setting values to NA/NULL to later be detected by the integrity check")
            parlance_name = ""
            species_prefix = ""
            taxon_classification = {}
            taxon_classification_check=False

        # Update species dictionary with new values
        logging.info("Updating keys for species table")
        species_dict['species'].update({
            'species_taxon_id': species_taxon_id,
            'parlance_name': parlance_name,
            'species_prefix': species_prefix})
        if taxon_classification_check:
            species_dict['taxonomy']= taxon_classification
            species_dict['taxonomy'].update({'lowest_taxon_id': species_dict['species']['lowest_taxon_id']})

        # Saving results
        output_file = os.path.basename(args.json_path).replace('.tmp', '.json')
        logging.info(f"Saving output: {output_file}")
        with open(output_file, 'w') as file:
            json.dump(species_dict, file)
        file.close()

    # Update taxonomy table, this is used as a standalone script
    if args.taxonomy_update and args.taxon_id:
        logging.info("Updating taxonomy table")
        taxon_dict = {}
        taxon_data = get_taxon_data(args.taxon_id , args.ncbi_url)
        taxon_classification = get_taxon_classification(taxon_data)
        taxon_dict['taxonomy']= taxon_classification
        taxon_dict['taxonomy'].update({'lowest_taxon_id': args.taxon_id})

        # Saving results
        output_file_taxon = f"taxonomy_{args.taxon_id}.json"
        logging.info(f"Saving output: {output_file_taxon}")
        with open(output_file_taxon, 'w') as file:
            json.dump(taxon_dict, file)
        file.close()

if __name__ == "__main__":
    main()