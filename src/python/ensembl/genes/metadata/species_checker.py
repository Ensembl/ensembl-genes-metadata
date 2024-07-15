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


import json
import argparse
import logging
import pymysql
import requests
import string
import random
import os

def species_taxon(lowest_taxon_id):
    
    uri = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/{lowest_taxon_id}/dataset_report"
    
    try:
        response = requests.get(uri)
        response.raise_for_status()
        
    except:
        response = requests.get(uri)
        response.raise_for_status()
        
    taxon_data = response.json()
    
    try:
        taxonomy = taxon_data['reports'][0]['taxonomy']['rank']
        if taxonomy in ['SUBSPECIES', 'STRAIN', 'VARIETAS', 'GENOTYPE', 'ISOLATE']:
            species_taxon_id = taxon_data['reports'][0]['taxonomy']['classification']['species']['id']
            logging.info(f"The assembly is a infraspecific taxon {taxonomy}")
        elif taxonomy == 'SPECIES':
            species_taxon_id = lowest_taxon_id
            logging.info(f"The assembly is a species taxon rank ")
        else:
            raise ValueError(f"Incorrect taxonomy ({taxonomy})")
    except KeyError:
        logging.info(f'Taxon do not have Rank available, retrieving information from another section of the report')
        species_taxon_id = taxon_data['reports'][0]['taxonomy']['classification']['species']['id']
    
    return species_taxon_id

def get_parlance_name(sci_name):
    
    parlance_file = "/Users/vianey/Documents/GitHub/core_meta_updates/scripts/metadata/snp_static.txt"
    data_dict = {}
    
    logging.info("Reading parlance name file (snp_static.txt) to look for a match")
    with open(parlance_file, 'r') as file:
        for line in file:
            key, value = line.rsplit("\t", 1)
            data_dict[key.strip()] = value.strip()
    
    parlance_name = data_dict.get(sci_name, "")
            
    return parlance_name

def getting_species_prefix(assembly_db):
    
    # Getting existing prefix from registry db
    conn = pymysql.connect(**assembly_db['registry'])
    cur  = conn.cursor()
    cur.execute("SELECT DISTINCT species_prefix FROM assembly")
    prefix_registry = cur.fetchall() 
    cur.close()
    
    # Getting existing prefix from metadata db
    conn = pymysql.connect(**assembly_db['metadata'])
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

def get_species_prefix(taxon_id, assembly_db):
    
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
        conn = pymysql.connect(**assembly_db['registry'])
        cur  = conn.cursor()
        query = f"SELECT DISTINCT species_prefix FROM assembly WHERE taxonomy = {taxon_id}"
        cur.execute(query)
        output_registry = cur.fetchall()
        cur.close()
        
        # Search prefix in DB (gb_assembly_metadata)
        conn = pymysql.connect(**assembly_db['metadata'])
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
            species_prefix = getting_species_prefix(assembly_db)
        # unique prefix detected
        elif len(prefix_list) == 1:
            species_prefix = prefix_list[0]
        # Multiple prefix detected, check species. 
        else:
            raise ValueError(f"The taxon {taxon_id} is already registered, multiple prefix detected: {prefix_list}")
        
    return species_prefix

def main():
    """Entry point
    """
    
    db_params = {
        "host":"mysql-ens-genebuild-prod-1",
        "user":"ensadmin",
        "password":"ensembl",
        "port": 4527,
        "database" : "gb_assembly_metadata_testing"}
    
    db_credentials = {
        'database': 'gb_assembly_registry',
        'host':'mysql-ens-genebuild-prod-1',
        'user':'ensro',
        'port': 4527 }
    
    assembly_db = {'registry': db_credentials, 'metadata': db_params}
    
    logging.basicConfig(filename="species_checker.log", level=logging.DEBUG, 
                        filemode='w', format="%(asctime)s:%(levelname)s:%(message)s")
    
    parser = argparse.ArgumentParser(prog="species_checker.py", 
                                    description="Update the keys of related to species table")
    
    parser.add_argument("--json-path", type=str,
                        help="Path to the Species's JSON to be updated")
    
    args = parser.parse_args()
    
    logging.info(f"Loading files: {args.json_path}")
    
    # Loading Species dictionary
    with open(args.json_path, 'r') as file:
        species_dict = json.load(file)
    file.close()
    
    # Retrieve missing information to add keys to species json
    logging.info("Updating keys for species table")
    species_taxon_id = species_taxon(species_dict['species']['lowest_taxon_id'])
    parlance_name = get_parlance_name(species_dict['species']['scientific_name'])
    species_prefix = get_species_prefix(species_dict['species']['lowest_taxon_id'], assembly_db)
    
    # Update species dictionary with new values
    species_dict['species'].update({
        'species_taxon_id': species_taxon_id,
        'parlance_name': parlance_name, 
        'species_prefix': species_prefix
        })
    
    # Saving results
    output_file = os.path.basename(args.json_path).replace('.tmp', '.json')
    logging.info(f"Saving output: {output_file}")
    with open(output_file, 'w') as file:
        json.dump(species_dict, file)
    file.close()

if __name__ == "__main__":
    main()