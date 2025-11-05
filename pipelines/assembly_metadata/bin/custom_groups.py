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

import argparse
import os
import pymysql
import json
import logging
from typing import Tuple

def connect_db(query:str, metadata_params) -> Tuple:
    logging.info(f"Querying metadata database with: {query}")
    conn = pymysql.connect(**metadata_params)
    with conn.cursor() as cursor:
        cursor.execute(query)
        result = cursor.fetchall()
    conn.close()
    return result

def get_group(asm_id:int, gca_chain:str, lowest_taxon:str, species_taxon:str, metadata_params) -> dict:
    """
    This function query the custom_group table to retrieve the group name for a given assembly or taxon.
    If no custom group available, it returns an empty dictionary.

    Args:
        asm_id (int): assembly_id of the assembly. Retrieved from the assembly table 
        gca_chain (str): assembly chain only
        lowest_taxon (str): species or subspecies taxon id. The lowest available for the taxon
        species_taxon (str): species taxon id of the assembly

    Raises:
        ValueError: If multiple custom groups are detected for the same assembly

    Returns:
        dict: Dictionary containing the group name and assembly id. Empty if no custom group is detected
    """
    
    logging.info(f"Retrieving custom group by assembly: {gca_chain}")
    asm_query = f"SELECT group_name FROM custom_group where item = '{gca_chain}'"
    rg1 = connect_db(asm_query, metadata_params)
    
    logging.info(f"Retrieving custom group by taxon: {species_taxon}")
    species_query = f"SELECT group_name FROM custom_group where item = '{species_taxon}'"
    rg2 =connect_db(species_query, metadata_params)

    taxon_query = f"SELECT group_name FROM custom_group where item = '{lowest_taxon}'"
    rg3 = connect_db(taxon_query, metadata_params)
    
    results = set([rg1, rg2, rg3])
    non_empty_results = [res for res in results if res]

    if len(non_empty_results) == 1:
        group_name = non_empty_results[0][0][0]
        logging.info("Custom group detected:", group_name)
        group_asm = {'group_assembly':
            {'assembly_id': asm_id, 'group_name': group_name}}
        
    elif len(non_empty_results) == 0:
        group_asm = {}
        logging.info("No custom group detected. Creating empty JSON file")
    else:
        logging.info("Multiple custom groups detected for the same assembly.")
        raise ValueError("Please check the custom group table")
    
    return group_asm   

def update_db(group_df: Tuple, metadata_params ) -> None:
    
    
    list_group_asm = []

    for group_name, group_type, item in group_df:
        # Assembly based groups
        if group_type == 'assembly':
            query = f"SELECT assembly_id FROM assembly where gca_chain = '{item}'"
            results = connect_db(query, metadata_params)
            
        # Taxon based groups
        if group_type == 'taxon':
            query = f"""
            SELECT assembly_id FROM assembly
            INNER JOIN species
            ON assembly.lowest_taxon_id = species.lowest_taxon_id
            WHERE (species.lowest_taxon_id = '{item}' OR species.species_taxon_id = '{item}')
            """
            results = connect_db(query, metadata_params)
        
        # Processing results, and filtering out the ones that are already in the group_assembly table
        for asm in results:
            assembly_id = asm[0]
            
            query = f"SELECT * FROM group_assembly WHERE group_name = '{group_name}' AND assembly_id = {assembly_id}"
            result= connect_db(query, metadata_params)

            if len(result) == 0:
                list_group_asm.append(
                            {'assembly_id': assembly_id,
                            'group_name': group_name})
                
    group_asm = {"group_assembly": list_group_asm}
        
    return group_asm


def main():
    """Module's entry point
    """

    logging.basicConfig(filename="custom_group.log", level=logging.DEBUG, filemode='w',
                        format="%(asctime)s:%(levelname)s:%(message)s")
    
    parser = argparse.ArgumentParser(prog='custom_group.py',
        description='Assign a custom group to an assembly or taxon when available')
    parser.add_argument('--accession', 
                        type=str,
                        help='Full GCA assembly accession')
    parser.add_argument('--update_db',
                        action='store_true',
                        help='Update previous loaded assemblies' )
    parser.add_argument('--metadata',
                        help = 'JSON file with metadata database parameters')
    
    args = parser.parse_args()
    logging.info(args)
    
    if os.path.exists(args.metadata):
        with open(args.metadata, "r") as file:
            metadata_params = json.load(file)
    else:
        raise FileNotFoundError("Metadata file not found")
    
    if args.accession and args.update_db:
        logging.info("Both accession and update_db flags are provided. Please provide only one")
        raise ValueError("Please provide only one flag")
    
    if args.accession: 
        logging.info(f"Accession provided: {args.accession}. Looking for match in custom_group table")
        
        gca_chain = args.accession.split('.')[0]    
        logging.info(f"Connecting to metadata database to retrieve data for {args.accession}")
        
        query = f"""
        SELECT asm.assembly_id, asm.lowest_taxon_id, asm.gca_chain, sp.species_taxon_id 
        FROM assembly asm
        INNER JOIN species sp
        on asm.lowest_taxon_id = sp.lowest_taxon_id 
        where asm.gca_chain = '{gca_chain}' ;"""
        asm_id, lowest_taxon, gca_chain, species_taxon = connect_db(query, metadata_params)[0]
        
        logging.info(f"The lowest taxon id is {lowest_taxon}, the species taxon is {species_taxon} and the gca chain is {gca_chain}")
        
        group_asm = get_group(asm_id, gca_chain, lowest_taxon, species_taxon, metadata_params)
        
        logging.info(f"Saving custom group json for {args.accession}")
        with open(f"{args.accession}_group.json", "w") as file:
            json.dump(group_asm, file)
        file.close()
        
    if args.update_db:
        # Querying the custom_group table
        item_query = f""" SELECT group_name, group_type, item FROM custom_group"""
        group_df = connect_db(item_query, metadata_params)
        
        # Getting values to update the group_assembly table
        group_asm = update_db(group_df, metadata_params)
        
        logging.info(f"Saving custom group json for general update")
        with open("update_db_group.json", "w") as file:
            json.dump(group_asm, file)
        file.close()
        
if __name__ == '__main__':
    main()    