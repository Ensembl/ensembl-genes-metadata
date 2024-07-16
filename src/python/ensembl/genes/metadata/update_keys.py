
import json
import os
import logging
import argparse

def update_keys(data_db, last_id_dict, table_conf):
    """Add keys to an already exiting dictionary

    Args:
        data_db (dict): It contains the data that intents to be insert/update to the DB
        last_id_dict (dict): It contain the table name and the respective row id value (primary key or foreign key)
        table_conf (dict): describe population method and dependant keys of the DB's table

    Returns:
        file: an updated json file containing the required keys
    """
    
    for table_name in data_db:
        logging.info(f"Table found in JSON: {table_name}")
        
        dkey = table_conf[table_name]['dkey']
        logging.info(f"Dependant key required for this table: {dkey}")
        
        dkey_value = last_id_dict[dkey]
        logging.info(f"Dependant key value: {dkey_value}")
        
        data_db[table_name].update({dkey:dkey_value})
        
    return data_db

def main():
    """Entry point"""
    
    table_conf = {
        "assembly": {"method": "per_col", "dkey":"None", "ukey": "None"},
        "organism": {"method": "per_col", "dkey":"None", "ukey": "None"},
        "species": {"method": "per_col", "dkey":"None", "ukey": "None"},
        "assembly_metrics" :  {"method": "per_row", "dkey":"assembly_id", "ukey": "None"},
        "bioproject_lineage": {"method": "per_row_key", "dkey":"assembly_id", "ukey": "None"},
        }
    
    logging.basicConfig(filename="update_keys.log", level=logging.DEBUG, filemode='w',
                    format="%(asctime)s:%(levelname)s:%(message)s")
    
    parser = argparse.ArgumentParser(prog="update_keys.py", 
                                    description="Update the keys of an existing json file")
    
    parser.add_argument("--json-path", type=str,
                        help="Path to the JSON to be updated")
    
    parser.add_argument("--file-id-path", type=str,
                        help="Path to the JSON file containing the keys to be added")
        
    # Parsing arguments 
    args = parser.parse_args()
    
    # Loading files
    logging.info(f"Loading files: {args.json_path}, {args.file_id_path}")
    
    with open(args.json_path, 'r') as file:
        data_db = json.load(file)
        file.close()

    with open(args.file_id_path, 'r') as file:
        last_id_dict = json.load(file)
        file.close()
    
    logging.info("Updating keys")
    data_db = update_keys(data_db, last_id_dict, table_conf)
    
    output_file = os.path.basename(args.json_path).replace('.tmp', '.json')
    
    logging.info(f"Saving output: {output_file}")
    with open(output_file, 'w') as file:
        json.dump(data_db, file)
    file.close()
    
if __name__ == '__main__':
    main()