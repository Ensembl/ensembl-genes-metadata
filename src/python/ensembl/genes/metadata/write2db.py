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
Script to create mysql queries (insert and update) given a json file as input
    
Raises:
    ValueError: rise and exception when key are missing of the data is not properly formatted

    Returns:
    str: MySQl queries
"""

import json
import argparse
import logging
import pymysql
import os

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

def insert_query(data_dict, table_name, table_conf, db_params):
    """This functions create a mysql insert queries using the input data provided

    Args:
        data_dict (dic): dictionary containing the key:values to be inserted
        table_name (_type_): table that is used to insert new data

    Returns:
        str: returns mysql query
    """
    
    if table_conf[table_name]['method'] in ['per_row', 'per_row_key']:
        logging.info(f"{table_name} is an attribute table (key:value pairs) ")
        
        dkey = table_conf[table_name]['dkey']
        
        # Getting columns names
        conn = pymysql.connect(**db_params)
        cur  = conn.cursor()
        cur.execute(f"SHOW COLUMNS FROM {table_name}")
        table_columns = cur.fetchall()
        columns = [column[0] for column in table_columns if column[3] != 'PRI']
        columns_string = ', '.join(columns)
        #
        value_list = []
        dkey_value = data_dict[dkey]
        
        for key,value in data_dict.items():
            if key !=  dkey: 
                if table_conf[table_name]['method'] in ['per_row']:
                    value_item = f"('{dkey_value}', '{key}', '{value}')"
                elif table_conf[table_name]['method'] in ['per_row_key']:
                    logging.info(f"{table_name} is an attribute table (key only) ")
                    value_item = f"('{dkey_value}', '{key}')"
                else:
                    raise ValueError(f"Invalid value in table config - method: {table_conf[table_name]['method'] } ")   
                
                value_list.append(value_item)
                values_string =  ', '.join(value_list)
        return f"""INSERT INTO {table_name} ({columns_string}) VALUES {values_string}"""

    else:
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
    
    return f"UPDATE {table_name} SET {update_values} WHERE {condition} ;"

def create_query(data_dict, table_name, update, table_conf, db_params):
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
        query = insert_query(data_dict, table_name, table_conf, db_params)
    
    return query

def retrieve_row_id(data_dict, table_name, db_params, table_conf):
    # Establishing connection to DB
    conn = pymysql.connect(**db_params)
    cur  = conn.cursor()

    # Getting id key name
    query = f"SHOW KEYS FROM {table_name} WHERE Key_name = 'PRIMARY'"
    cur.execute(query)
    id_name = cur.fetchone()[4]
    logging.info(f"Retrieving value of {id_name} from {table_name}")

    # Getting constraint keys
    constraint_query = f"""SELECT column_name FROM information_schema.key_column_usage
    WHERE
        table_schema = 'gb_assembly_metadata_testing'
        AND table_name = '{table_name}'
        AND constraint_name != 'PRIMARY'
        AND referenced_table_name IS NULL;"""
    cur.execute(constraint_query)
    constraint= cur.fetchall()
    logging.info(f" Detected uniqueness constrains: {constraint}")
    
    # Per row method
    if table_conf[table_name]['method'] in ['per_row', 'per_row_key']:
        print(f"{table_name} is an attribute table (key:value paris) ")
        dkey = table_conf[table_name]['dkey']
        
        # Getting columns names
        conn = pymysql.connect(**db_params)
        cur  = conn.cursor()
        cur.execute(f"SHOW COLUMNS FROM {table_name}")
        table_columns = cur.fetchall()
        columns = [column[0] for column in table_columns if column[3] != 'PRI']
        
        dkey_value = data_dict[dkey]
        
        for key,value in data_dict.items():
            if key !=  dkey: 
                
                if table_conf[table_name]['method'] in ['per_row']:
                    condition_string = f"{columns[0]} = '{dkey_value}' AND {columns[1]} =  '{key}' AND {columns[2]} = '{value}'"
                elif table_conf[table_name]['method'] in ['per_row_key']:
                    condition_string = f"{columns[0]} = '{dkey_value}' AND {columns[1]} =  '{key}'"
                else:
                    raise ValueError(f"Invalid value in table config - method: {table_conf[table_name]['method'] } ")
                
                conn = pymysql.connect(**db_params)
                cur  = conn.cursor()
                retrieving_query = f"SELECT {id_name} FROM {table_name} WHERE {condition_string} ;" 
                cur.execute(retrieving_query)
                last_id_tmp = cur.fetchall()
                cur.close()

                if len(last_id_tmp) > 1:
                    raise ValueError(f"The query retrieves more than more value, unique value expected {retrieving_query}")
                else:
                    last_id = last_id_tmp[0][0]
    
    elif table_conf[table_name]['method'] == 'per_col':
        # Building conditionals based on uniqueness constrain 
        condition_list = []
        for key in constraint:
            condition_list.append(f"{key[0]} = '{data_dict[key[0]]}'")

        condition_string = ' AND '.join(condition_list)

        retrieving_query = f"SELECT {id_name} FROM {table_name} WHERE {condition_string} ;" 
        cur.execute(retrieving_query)
        last_id_tmp = cur.fetchall()
        last_id_tmp

        if len(last_id_tmp) > 1:
            raise ValueError(f"The query retrieves more than more value, unique value expected {retrieving_query}")
        else:
            last_id = last_id_tmp[0][0]
        
    else:
        raise ValueError(f"Failed check of {table_name} configuration")
        
    return last_id

def execute_query(query, db_params, table_name, data_dict, table_conf):
    # Connecting to db 
    conn = pymysql.connect(**db_params)
    cur  = conn.cursor()
    # Getting id name
    cur.execute(f"SHOW KEYS FROM {table_name} WHERE Key_name = 'PRIMARY'")
    id_name = cur.fetchone()[4]
    
    try:
        # Execute query 
        cur.execute(query)
        # Getting id 
        id_value = cur.lastrowid
        
    except pymysql.IntegrityError as e:
        if e.args[0] == 1062:
            # Retrieve value of a already exiting row
            id_value = retrieve_row_id(data_dict, table_name, db_params, table_conf)
        
    except Exception as ee:
        print(f'Error: {ee}')
        raise ValueError
    
    cur.close()
    conn.close()
    
    return id_value, id_name

def main():
    """Entry point"""
    
    table_conf = {
        "assembly": {"method": "per_col", "dkey":"None", "ukey": "None"},
        "organism": {"method": "per_col", "dkey":"None", "ukey": "None"},
        "species": {"method": "per_col", "dkey":"None", "ukey": "None"},
        "assembly_metrics" :  {"method": "per_row", "dkey":"assembly_id", "ukey": "None"},
        "bioproject_lineage" :  {"method": "per_row_key", "dkey":"assembly_id", "ukey": "None"},
        }
    
    db_params = {
        "host":"mysql-ens-genebuild-prod-1",
        "user":"ensadmin",
        "password":"ensembl",
        "port": 4527,
        "database" : "gb_assembly_metadata_testing"
    }
    
    logging.basicConfig(filename="write2db.log", level=logging.DEBUG, filemode='w',
                    format="%(asctime)s:%(levelname)s:%(message)s")
    
    parser = argparse.ArgumentParser(prog="create_query.py", 
                                    description="Create insert/update queries and execute them in the target DB")
    
    parser.add_argument("--file-path", type=str,
                        help="Path to the JSON file containing data to insert/update in a DB")
    
    #parser.add_argument('--table-config', type=str,
    #                    help="file with the configuration of the DB")
    
    #parser.add_argument('--db', type=str,
    #                    help="file with the credentials of the DB")
    
    parser.add_argument('--update', action='store_true',
                        help="If option is added it indicates the input data will be used to update")
    
    # Parsing arguments 
    args = parser.parse_args()
    
    # Loading files
    logging.info(f"Loading file: {args.file_path}")
    file = open(str(args.file_path))
    input_data = json.load(file)
    file.close()
    
    # Update optional command 
    update = args.update
    
    # Output
    root_name, _ = os.path.splitext(args.file_path) 
    output = root_name + ".last_id"
    with open(output, 'w') as file:
        pass
    
    last_id_dict = {}
    
    # Module 
    for table_name in input_data:
        logging.info(f"Processing input data, loading {table_name} table")
        # Check input data structure
        check = check_dict_structure(input_data[table_name])
        #logging.info(check)
        if check:
            logging.info("Lists of dictionaries detected, processing each dictionary")
            for row in input_data[table_name]:
                query = create_query(row, table_name, update, table_conf, db_params)
                id_value, id_name= execute_query(query,db_params, table_name, input_data[table_name], table_conf)
                logging.info(f"Data was inserted in {table_name}. Last value of {id_name} is {id_value}")
                # saving last id in dict
                last_id_dict.update({id_name:id_value})
                
        else:
            logging.info("Regular dictionary detected, processing key:value pair values")
            # Data is a dictionary (This part is not tested yet)
            query = create_query(input_data[table_name], table_name, update, table_conf, db_params)
            id_value, id_name= execute_query(query,db_params,table_name, input_data[table_name], table_conf)
            logging.info(f"Data was inserted in {table_name}. Last value of {id_name} is {id_value}")
            # saving last id in dict
            last_id_dict.update({id_name:id_value})
    
    # saving output       
    with open(output, 'a') as file:
        json.dump(last_id_dict, file)
    file.close()

if __name__ == "__main__":
    main()