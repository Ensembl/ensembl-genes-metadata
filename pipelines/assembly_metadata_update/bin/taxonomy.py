#!/usr/bin/env python3

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
import json
import logging
from typing import Dict
import pymysql # type: ignore


def execute_query(query, db_params):
    conn = pymysql.connect(**db_params)
    cursor = conn.cursor()
    cursor.execute(query)
    result = cursor.fetchall()
    cursor.close()
    conn.close()
    return result

def comparing_basic_taxon_data(data, accession, metadata_params):

    output_line_list = []

    # Getting info from NCBI report
    taxon_id_ncbi = data['reports'][0].get('organism', '').get('tax_id')
    organism_name_ncbi = data['reports'][0].get('organism', '').get('organism_name')
    common_name_ncbi = data['reports'][0].get('organism', '').get('common_name')
    

    # Taxon ID check
    query_taxon_id = f"SELECT lowest_taxon_id from assembly WHERE CONCAT(gca_chain, '.', gca_version) = '{accession}'"
    taxon_id = execute_query(query_taxon_id, metadata_params)[0][0]

    if taxon_id_ncbi == taxon_id:
        logging.info(f"No update needed for taxon ID of assembly {accession}")
    else:
        logging.info(f"Update needed for taxon ID of assembly {accession}: current taxon ID in Registry is {taxon_id}, taxon ID from NCBI is {taxon_id_ncbi}")
        query_update_taxon_id = f"""UPDATE assembly 
        SET lowest_taxon_id = '{taxon_id_ncbi}' 
        WHERE CONCAT(gca_chain, '.', gca_version) = '{accession}';"""
        logging.info(query_update_taxon_id)
        output_line = f"{accession}, taxon_id_check, true, {taxon_id}, {taxon_id_ncbi}"
        output_line_list.append(output_line)


    # Scientific and common name check
    query_species = f"SELECT scientific_name, common_name from species where lowest_taxon_id = '{taxon_id_ncbi}';"
    scientific_name, common_name = execute_query(query_species, metadata_params)[0]

    if scientific_name != organism_name_ncbi:
        logging.info(f"Update required for scientific name of taxon ID {taxon_id_ncbi}: current scientific name in Registry is {scientific_name}, scientific name from NCBI is {organism_name_ncbi}")
        query_update_scientific_name = f"""UPDATE species 
        SET scientific_name = '{organism_name_ncbi}' 
        WHERE lowest_taxon_id = '{taxon_id_ncbi}';"""
        logging.info(query_update_scientific_name)
        output_line = f"{accession}, scientific_name_check, true, {scientific_name}, {organism_name_ncbi}"
        output_line_list.append(output_line)
    else:
        logging.info(f"No update required for scientific name of taxon ID {taxon_id_ncbi}")
    
    if common_name != common_name_ncbi:
        logging.info(f"Update required for common name of taxon ID {taxon_id_ncbi}: current common name in Registry is {common_name}, common name from NCBI is {common_name_ncbi}")
        query_update_common_name = f"""UPDATE species 
        SET common_name = '{common_name_ncbi}' 
        WHERE lowest_taxon_id = '{taxon_id_ncbi}';"""
        logging.info(query_update_common_name)
        output_line = f"{accession}, common_name_check, true, {common_name}, {common_name_ncbi}"
        output_line_list.append(output_line)
    else:
        logging.info(f"No update required for common name of taxon ID {taxon_id_ncbi}")

    return output_line_list
    
def main():
    """ Module's entry point
    """

    logging.basicConfig(filename="update_taxonomy.log", level=logging.DEBUG, filemode='w',
                    format="%(asctime)s:%(levelname)s:%(message)s")

    parser = argparse.ArgumentParser(prog='taxonomy.py',
                                    description="Retrieve metadata from NCBI API for a given GCA accession and store it in JSON files to be inserted in the database.")

    parser.add_argument('--accession_json',
                        type=str,
                        required=True,
                        help='Path to JSON file with assembly metadata retrieved from NCBI API')
    parser.add_argument('--accession',
                        type=str,
                        required=True,
                        help='GCA accession to retrieve metadata')
    parser.add_argument('--metadata_params',
                        type=str,
                        required=True,
                        help='Database connection parameters for metadata database in JSON format')
    

    args = parser.parse_args()
    logging.info(args)

    accession = args.accession.strip()

    with open(args.accession_json, 'r') as json_file:
        data = json.load(json_file)
    
    with open(args.metadata_params, 'r') as params_file:
        metadata_params = json.load(params_file)

    output_line_list = comparing_basic_taxon_data(data, accession, metadata_params)

    for output_line in output_line_list:
        print(output_line)

if __name__ == '__main__':
    main()
