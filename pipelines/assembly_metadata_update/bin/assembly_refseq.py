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

def comparing_refseq(data, accession, metadata_params):
    """Compare available paired REFSEQ from NCBI and Registry and generate update queries if needed.
    Args:
        data (dict): JSON data retrieved from NCBI API for the assembly
        accession (str): GCA accession of the assembly
        metadata_params (dict): Database connection parameters for metadata database
    Returns:
        str: A line with accession, current refseq in Registry, and refseq from NCBI
    """

    # Getting info from Registry
    query_refseq = f"SELECT refseq_accession FROM assembly WHERE CONCAT(gca_chain, '.', gca_version) = '{accession}'"
    refseq_accession = execute_query(query_refseq, metadata_params)[0][0]

    # Getting info from NCBI
    paired_accession =  data['reports'][0].get('paired_accession',"")

    # Comparison
    if refseq_accession == paired_accession:
        logging.info(f"No update needed for RefSeq accession of assembly {accession}")
        output_line = f"{accession}, refseq_check, false, NA, NA"
    elif not refseq_accession and paired_accession != "":
        logging.info(f"No RefSeq accession found for assembly {accession} in Registry, setting to {paired_accession}")
        query_update_refseq = f"UPDATE assembly a SET refseq_accession = '{paired_accession}' WHERE CONCAT(a.gca_chain, '.', a.gca_version) = '{accession}';"
        logging.info(query_update_refseq)
        affected = execute_write(query_update_refseq, metadata_params)
        output_line = f"{accession}, refseq_check, true, no_refseq, {paired_accession}"
    elif not refseq_accession and paired_accession == "":
        logging.info(f"No update needed for RefSeq accession of assembly {accession}")
        output_line = f"{accession}, refseq_check, false, NA, NA"
    else:   
        logging.info(f"Update needed for RefSeq accession of assembly {accession}: current RefSeq in Registry is {refseq_accession}, RefSeq from NCBI is {paired_accession}")
        query_update_refseq = f"UPDATE assembly a SET refseq_accession = '{paired_accession}' WHERE CONCAT(a.gca_chain, '.', a.gca_version) = '{accession}';"
        logging.info(query_update_refseq)
        affected = execute_write(query_update_refseq, metadata_params)
        output_line = f"{accession}, refseq_check, true, {refseq_accession}, {paired_accession}"

    return output_line

    
def main():
    """ Module's entry point
    """

    logging.basicConfig(filename="update_assembly_refseq.log", level=logging.DEBUG, filemode='w',
                    format="%(asctime)s:%(levelname)s:%(message)s")

    parser = argparse.ArgumentParser(prog='update_assembly_refseq.py',
                                    description="Checks if there is an new available paired refseq accession.")

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

    output_line = comparing_refseq(data, accession, metadata_params)
    print(output_line)

if __name__ == '__main__':
    main()
