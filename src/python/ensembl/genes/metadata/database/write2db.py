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


#import pymysql
import json
import argparse
import logging

def check_dict_structure(input_dict):
    """This functions checks the structure of the database, 
    identify if the data is a dictionary or a list of dictionary
    
    Args:
        input_dict (dict or list[dict]): input data to load in db

    Returns:
        boolean: returns a True if the input data is a list of dictionaries
                    or False if the input data is a dictionary
    """
    
    if isinstance(input_dict, list):
        dict_list = True
    else:
        for key, value in input_dict.items():
            if isinstance(value, list):  # Check if the value is a list
                if all(isinstance(item, dict) for item in value):  # Check if all items in the list are dictionaries
                    #print(f"Key '{key}' is linked to a list of dictionaries.")
                    dict_list = True
                else:
                    raise ValueError(f"Table '{key}' is not properly formatted")
            else:
                #print(f"Key '{key}' is not linked to a list of dictionaries.")
                dict_list = False
            
    return dict_list

def check_key(data_dict, table_name, update, table_conf):
    """
    This function check if the data have all keys necessary for a successful execution.
    It is used for insert queries or update queries

    Args:
        data_dict (dict): input data of one table to load in db
        table_name (str): table to load in db
        update (boolean): True if the data is to update an exiting row. False is default (insert query)

    Raises:
        ValueError: raise an error when keys are missing
    """
    
    if update: #check update key -> ukey
        key = 'ukey'
    else: # check dependent key -> dkey
        key = 'dkey'
        
    # If insert query do not require dkey:
    if table_conf[table_name][key] == 'None':
        logging.info(f"The {table_name} table does not require any dependent/update key. Update {update} ")
        
    # If insert query do require dkey:
    elif table_conf[table_name][key] != 'None':
        if table_conf[table_name][key] in data_dict.keys():
            logging.info(f'Key {table_conf[table_name][key]} exist in {table_name}. Update {update}')
        else:
            raise ValueError(f"key {table_conf[table_name][key]} not found in {table_name} table. Update {update}")
    else:
        raise ValueError(f"Unexpected value {table_conf[table_name][key]} for table {table_name} table. Update {update}")

def insert_query(data_dict, table_name):
    """This functions create a mysql insert queries using the input data provided

    Args:
        data_dict (dic): dictionary containing the key:values to be inserted
        table_name (_type_): table that is used to insert new data

    Returns:
        str: returns mysql query
    """
    # crete basic query
    table_var_string = ", ".join(list(data_dict.keys()))
    values_strings = ','.join([f"'{value}'" for value in list(data_dict.values())]).replace("''" , "NULL" )
    
    return f"""INSERT INTO {table_name} ({table_var_string}) VALUES ({values_strings}) ;"""

def update_query(data_dict, table_name, table_conf):
    """
    This functions create a mysql update query using the input data provided

    Args:
        data_dict (dict): dictionary containing the key:values to be inserted
        table_name (str): table that is being updated

    Returns:
        str: returns a mysql query 
    """
    update_list = []
    for key, value in data_dict.items():
        if table_conf[table_name]['ukey'] != key:
            update_list.append(f"{key} = {value}")
        else:
            condition = f"{key} = {value}"
        
    update_values = ','.join(update_list)
    
    return f"UPDATE {table_name} SET {update_values} WHERE {condition}"

def create_query(data_dict, table_name, update, table_conf):
    """
    This function create an insert or update MySQL query depending if update argument was provided (True).
    It use the check_key function to determinate of the data provided has the enough keys to create the query.
    
    Args:
        data_dict (dict): dictionary containing key:values to be insert/update in the DB
        table_name (str): table name to be used for the insert/update operation
        update (bool): boolean variable indicating if the operation is an update

    Returns:
        str: mysql query 
    """
    # checking if relevant keys are missing
    check_key(data_dict, table_name, update, table_conf)
    
    if update: #input data will be used to update a row
        query = update_query(data_dict, table_name, table_conf)
    else: # input data will be used to insert a new row
        query = insert_query(data_dict, table_name)
    
    return query

def main():
    """Entry point"""
    
    table_conf = {"run": {"method": "per_col", "dkey":"None", "ukey": "None"},
            "study": {"method": "per_col", "dkey":"None", "ukey": "None"},
            "data_files": {"method": "per_col","dkey": "run_id", "ukey": "file_id"},
            "align": {"method": "per_col",  "dkey": "run_id", "ukey": "None"}
        }

    db_params = {
        "host":"mysql-ens-genebuild-prod-1",
        "user":"ensadmin",
        "password":"ensembl",
        "port": "4527",
        "database" : "transcriptomic_hackathon_new"
    }
    
    logging.basicConfig(filename="write2db.log", level=logging.DEBUG, filemode='w',
                    format="%(asctime)s:%(levelname)s:%(message)s")
    
    parser = argparse.ArgumentParser(prog="write.py", 
                                    description="Write JSON files to DB")
    
    parser.add_argument("--file-path", type=str,
                        help="JSON file containing data to insert/update in a DB")
    
    #parser.add_argument('--table-config', type=str,
    #                    help="file with the configuration of the DB")
    
    #parser.add_argument('--db', type=str,
    #                    help="file with the credentials of the DB")
    
    parser.add_argument('--update', action='store_true',
                        help="Set to True if the input data will be used to update")
    
    # Parsing arguments 
    args = parser.parse_args()
    #print(args)
    file = open(str(args.file_path))
    input_data = json.load(file)
    file.close()
    update = args.update
    
    # Module 
    for table_name in input_data:
        logging.info(f"Processing input data, loading {table_name} table")
        # Check input data structure
        check = check_dict_structure(input_data[table_name])
        #print(check)
        if check:
            logging.info("List of dictionaries detected, processing each dictionary")
            for row in input_data[table_name]:
                query = create_query(row, table_name, update, table_conf)
                print(query)
        else:
            # Data is a dictionary
            query = create_query(input_data[table_name], table_name, update, table_conf)
            print(query)

if __name__ == "__main__":
    main()
