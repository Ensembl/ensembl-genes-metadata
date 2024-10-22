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
import pymysql
import logging
import json
import os

def delete_assembly(assembly_id: int, metadata_params: dict):
    """
    Delete records from metadata database for a given assembly_id.
    Affects assembly, assembly_metrics, organism and bioproject tables

    Args:
        assembly_id (int): Assembly ID to be deleted
    """
    con = pymysql.connect(**metadata_params)  
    cur = con.cursor()
    
    # Delete records from all metadata tables
    cur.execute(f"DELETE FROM assembly WHERE assembly_id = '{assembly_id}'")
    cur.execute(f"DELETE FROM assembly_metrics WHERE assembly_id = '{assembly_id}'")
    cur.execute(f"DELETE FROM organism WHERE assembly_id = '{assembly_id}'")
    cur.execute(f"DELETE FROM bioproject WHERE assembly_id = '{assembly_id}'")
    print(f"Records for assembly_id {assembly_id} were deleted")
    con.close()

def integrity_mode(delete, metadata_params):
    
    logging.info("Connecting to Assembly metadata database to retrieve a list of assembly_ids ")
    con = pymysql.connect(**metadata_params)  
    cur = con.cursor()
    cur.execute("SELECT assembly_id, CONCAT(GCA_chain, '.', gca_version) FROM assembly")
    gca_list = cur.fetchall()
    con.close()
    print(f"Number of assemblies in the database: {len(gca_list)}")
    
    for assembly_id, gca in gca_list:
        #logging.info(f"Checking records for assembly_id {assembly_id} in assembly_metrics, organism and bioproject tables")
        con = pymysql.connect(**metadata_params)  
        cur = con.cursor()
        
        cur.execute(f"SELECT * FROM assembly_metrics WHERE assembly_id = '{assembly_id}'")
        data = cur.fetchall()
        metrics_table = len(data) > 1
        
        cur.execute(f"SELECT * FROM organism WHERE assembly_id = '{assembly_id}'")
        data = cur.fetchall()
        organism_table = len(data) == 1
        
        cur.execute(f"SELECT * FROM bioproject WHERE assembly_id = '{assembly_id}'")
        data = cur.fetchall()
        bioproject_table = len(data) >= 1
        
        con.close()
        
        if not (metrics_table and organism_table and bioproject_table):
            msg = f"Accession {gca} has missing data. Remove records associate to assembly_id {assembly_id}"
            print(msg)
            logging.info(msg)
            if delete:
                delete_assembly(assembly_id, metadata_params)

def per_gca_mode(gca_list, metadata_params):
    for gca in gca_list:
        logging.info(f"Connecting to metadata database to retrieve assembly_id for {gca}")
        con = pymysql.connect(**metadata_params)  
        cur = con.cursor()
        query = f"SELECT assembly_id FROM assembly WHERE CONCAT(GCA_chain, '.', gca_version) = '{gca}'"
        #logging.info(query)
        cur.execute(query)
        result = cur.fetchall()
        con.close()
        
        #logging.info(f"Results from the query: {result}")
        if len(result) != 0:
            # First check
            if len(result) > 1:
                raise ValueError("There are more than one assembly_id associated to this GCA")
            
            # Second check
            if len(result[0]) != 1:
                raise ValueError(f"Incorrect data retrieved for {gca}. Only assembly id, tuple should be length 1")
            
            # If all checks are passed, delete records
            asm_id = result[0][0]
            logging.info(f"Deleting records for assembly_id {asm_id}: {gca}")
            delete_assembly(asm_id, metadata_params)
            
        else:
            msg = f"No records found for {gca} in assembly table"
            logging.info(msg)
            print(msg)
            pass

def main():
    """
    Module's entry point
    """
    
    logging.basicConfig(filename="clean_gca_records.log", level=logging.DEBUG, filemode='w',
                        format="%(asctime)s:%(levelname)s:%(message)s")
    
    parser = argparse.ArgumentParser(prog='clean_gca_records.py',
        description='Clean records from metadata database based on GCA records completeness o a list of GCA accessions')
    parser.add_argument('--file-list', 
                        help='TXT file with a list of GCA accessions')
    parser.add_argument('--integrity', 
                        action='store_true',
                        help='Check integrity of the records. Reports records with missing data')
    parser.add_argument('--delete', 
                        default=False,
                        action='store_true',
                        help='Delete records with missing data')
    parser.add_argument('--metadata',
                        help='Path to metadata database connection parameters')
    
    args = parser.parse_args()
    logging.info(args)
    
    if args.metadata:
        if not os.path.exists(args.metadata):
            raise ValueError("Metadata params json file does not exist")
        else:
            with open(args.metadata, 'r') as f:
                metadata_params = json.load(f)
                f.close()
    
    if not args.integrity and not args.file_list:
        raise ValueError("No action was provided")
    
    if args.integrity:
        logging.info("Checking integrity of the records")
        integrity_mode(args.delete, metadata_params)
    
    if args.file_list:
        logging.info(f"Reading GCA list from {args.file_list}")
        with open(args.file_list, 'r') as f:
            gca_list = f.read().splitlines()
            logging.info(gca_list)
            f.close()
        
        msg = f"Number of GCAs in the list: {len(gca_list)}"
        logging.info(msg)
        print(msg)
        per_gca_mode(gca_list, metadata_params)

if __name__ == "__main__":
    main()