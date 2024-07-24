
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

"""The module will update the keys of a .tmp file. The keys are extracted from a JSON file created by write2db.py. 
    The new JSON file are saved with the same name as the input file but with the .json extension.
    
Args:
    json_path (str): Path to the JSON-like file with .tmp extension to be updated 
    file_id_path (str): Path to the JSON file containing the keys to be added

Returns:
    str: a JSON file with the updated keys
"""

import json
import os
import logging
import argparse
from typing import Dict, Any
from db_table_conf import TABLE_CONF

def update_keys(data_db: Dict[str, Any], last_id_dict: Dict[str, Any]) -> Dict[str, Any]: 
    """Insert the missing keys in the data_db dictionary using the last_id_dict dictionary.

    Args:
        data_db (dict): data that intents to be insert/update to the DB and is missing a key
        last_id_dict (dict): dictionary with the table name as key and the row id value as value 

    Returns:
        dict: an updated dictionary containing the data and required keys
    """
    
    for table_name in data_db:
        logging.info(f"Table found in JSON: {table_name}")
        
        dkey = TABLE_CONF[table_name]['dkey']
        logging.info(f"Dependant key required for this table: {dkey}")
        
        dkey_value = last_id_dict[dkey]
        logging.info(f"Dependant key value: {dkey_value}")
        
        data_db[table_name].update({dkey:dkey_value})
        
    return data_db

def main():
    """Module's entry point"""
    
    logging.basicConfig(filename="update_keys.log", level=logging.DEBUG, filemode='w',
                    format="%(asctime)s:%(levelname)s:%(message)s")
    
    parser = argparse.ArgumentParser(prog="update_keys.py", 
                                    description="Update keys in a JSON-like file. The keys are extracted from a JSON file created by write2db.py. The new JSON file is saved with the same name as the input file but with the .json extension.")
    
    parser.add_argument("--json-path", type=str,
                        help="Path to the JSON-like file (.tmp) to be updated")
    
    parser.add_argument("--file-id-path", type=str,
                        help="Path to the JSON file containing the keys to be added. File is the output of write2db.py")
        
    # Parsing arguments 
    args = parser.parse_args()
    logging.info(args)
    
    # Loading files
    logging.info(f"Loading files: {args.json_path}, {args.file_id_path}")
    
    with open(args.json_path, 'r') as file:
        data_db = json.load(file)
        file.close()

    with open(args.file_id_path, 'r') as file:
        last_id_dict = json.load(file)
        file.close()
    
    logging.info("Updating keys")
    data_db = update_keys(data_db, last_id_dict)
    
    output_file = os.path.basename(args.json_path).replace('.tmp', '.json')
    
    logging.info(f"Saving output: {output_file}")
    with open(output_file, 'w') as file:
        json.dump(data_db, file)
    file.close()
    
if __name__ == '__main__':
    main()