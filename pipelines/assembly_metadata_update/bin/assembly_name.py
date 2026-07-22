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
from typing import Dict, Any
import pymysql


def execute_query(query, db_params):
    conn = pymysql.connect(**db_params)
    cursor = conn.cursor()
    cursor.execute(query)
    result = cursor.fetchall()
    cursor.close()
    conn.close()
    return result

def execute_write(query: str, db_params: Dict[str, Any]) -> int:
    """
    Execute INSERT/UPDATE/DELETE and commit.
    Returns number of affected rows.
    """
    conn = pymysql.connect(**db_params)
    try:
        with conn.cursor() as cursor:
            affected = cursor.execute(query)
        conn.commit()
        return affected
    finally:
        conn.close()

def comparing_name(data, accession, metadata_params):
    """Compare assembly status from NCBI and Registry and generate update queries if needed.
    Args:
        data (dict): JSON data retrieved from NCBI API for the assembly
        accession (str): GCA accession of the assembly
        metadata_params (dict): Database connection parameters for metadata database
    Returns:
        str: A line with accession, current status in Registry, and status from NCBI
    """

    # Getting info from Registry
    query_name = f"SELECT asm_name FROM assembly WHERE CONCAT(gca_chain, '.', gca_version) = '{accession}'"
    asm_name_registry = execute_query(query_name, metadata_params)[0][0]

    # Getting info from NCBI
    asm_name_ncbi = data['reports'][0].get('assembly_info').get('assembly_name')

    # Comparing and generating update queries
    if asm_name_registry == asm_name_ncbi:
        logging.info(f"No update needed for assembly {accession}. Current name in Registry and NCBI is {asm_name_registry}")
        output_line = f"{accession}, asm_name_check, false, NA, NA"
    else:
        logging.info(f"Update needed for assembly {accession}: current name in Registry is {asm_name_registry}, name from NCBI is {asm_name_ncbi}")
        query_update_name = f"UPDATE assembly SET asm_name = '{asm_name_ncbi}' WHERE CONCAT(gca_chain, '.', gca_version) = '{accession}';"
        logging.info(query_update_name)
        affected = execute_write(query_update_name, metadata_params)
        output_line = f"{accession}, asm_name_check, true, {asm_name_registry}, {asm_name_ncbi}"

    return output_line

    
def main():
    """ Module's entry point
    """

    logging.basicConfig(filename="update_assembly_name.log", level=logging.DEBUG, filemode='w',
                    format="%(asctime)s:%(levelname)s:%(message)s")

    parser = argparse.ArgumentParser(prog='assembly_name.py',
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

    output_line = comparing_name(data, accession, metadata_params)
    print(output_line)

if __name__ == '__main__':
    main()
