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
This module is used to create insert or update queries for a MySQL database. 
The module execute the queries in the target database and return the id of the row inserted or updated.

Args:
    file_path (str): Path to the JSON file containing data to insert or update in a DB

Raises:
    ValueError: when the query is not executed successfully
    ValueError: when keys are missing
    ValueError: when the input data is not properly formatted
    ValueError: when there are invalid method values in the table configuration
    ValueError: when the query retrieves more than one value
    ValueError: when table name is not found in configuration

Returns:
    str: json file with last id of the row inserted or updated
"""

import json
import argparse
import logging
import pymysql
import os
from typing import Dict, Tuple, Any

def check_dict_structure(input_dict) -> bool:
    """This functions checks the structure of the dictionary to,
    identify if the data is a dictionary or a list of dictionary

    Args:
        input_dict (dict or list[dict]): input data to load in db

    Returns:
        boolean: returns a True if the input data is a list of dictionaries
                    or False if the input data is a dictionary
    """
    
    if isinstance(input_dict, list):
        dict_islist = True
    else:
        for key, value in input_dict.items():
            if isinstance(value, list):  # Check if the value is a list
                if all(isinstance(item, dict) for item in value):  # Check if all items in the list are dictionaries
                    #print(f"Key '{key}' is linked to a list of dictionaries.")
                    dict_islist = True
                else:
                    raise ValueError(f"Table '{key}' is not properly formatted")
            else:
                #print(f"Key '{key}' is not linked to a list of dictionaries.")
                dict_islist = False
            
    return dict_islist

def check_key(data_dict, table_name, update, table_conf) -> None:
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

def insert_query(data_dict: Dict , table_name: str, table_conf, metadata_params) -> str:
    """This functions create a mysql insert queries using the input data provided

    Args:
        data_dict (dic): dictionary containing the key:values to be inserted
        table_name (str): table that is used to insert new data

    Returns:
        str: returns mysql query
    """
    
    if table_conf[table_name]['method'] in ['per_row', 'per_row_key']:
        logging.info(f"{table_name} is an attribute table (key:value pairs) ")
        
        dkey = table_conf[table_name]['dkey']
        
        # Getting columns names
        conn = pymysql.connect(**metadata_params)
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

def update_query(data_dict: Dict, table_name: str, table_conf) -> str:
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
        if table_conf[table_name]['ukey'] == key:
            condition = f"{key} = {value}"
        else:
            update_list.append(f"{key} = '{value}'")

    update_values = ','.join(update_list)

    return f"UPDATE {table_name} SET {update_values} WHERE {condition} ;"

def create_query(data_dict: Dict, table_name: str, update: bool, table_conf, metadata_params) -> str:
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
        query = insert_query(data_dict, table_name, table_conf, metadata_params)
    
    logging.info(f"Query created: {query}")
    return query

def retrieve_row_id(data_dict: Dict, table_name: str, table_conf, metadata_params) -> int:
    """
    If data is already inserted in the db, this function retrieves the id of the row using the data provided and the table name, 
    it will retrieve the uniqueness constrain of the table and use it to retrieve the id of the row.

    Args:
        data_dict (dict): dictionary containing key:values to be insert/update in the DB
        table_name (str): table name to be used for the insert/update operation

    Returns:
        int: id of the row
    """

    # Establishing connection to DB
    conn = pymysql.connect(**metadata_params)
    cur  = conn.cursor()

    # Getting id key name
    query = f"SHOW KEYS FROM {table_name} WHERE Key_name = 'PRIMARY'"
    cur.execute(query)
    id_name = cur.fetchone()[4]
    logging.info(f"Retrieving value of {id_name} from {table_name}")

    # Getting constraint keys
    constraint_query = f"""SELECT column_name FROM information_schema.key_column_usage
    WHERE
        table_schema = 'gb_assembly_metadata'
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
        conn = pymysql.connect(**metadata_params)
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

                conn = pymysql.connect(**metadata_params)
                cur  = conn.cursor()
                retrieving_query = f"SELECT {id_name} FROM {table_name} WHERE {condition_string} ;"
                logging.info(f"Retrieving query: {retrieving_query}")
                cur.execute(retrieving_query)
                last_id_tmp = cur.fetchall()
                cur.close()

                if len(last_id_tmp) > 1:
                    raise ValueError(f"The query retrieves more than more value, unique value expected {retrieving_query}")
                elif last_id_tmp == () and table_conf[table_name]['method'] in ['per_row']:
                    logging.info("Failed to retrieve value for last id. Inserting missing data")
                    query_missing_insert = f"INSERT INTO {table_name} ({columns[0]}, {columns[1]}, {columns[2]}) VALUES ('{dkey_value}', '{key}', '{value}') ;"
                    logging.info("Insert query: %s", query_missing_insert)
                    conn = pymysql.connect(**metadata_params)
                    cur  = conn.cursor()
                    cur.execute(query_missing_insert)
                    last_id = cur.lastrowid
                    conn.close()
                else:
                    logging.info(f"Retrieved value for last id: {last_id_tmp}")
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
        elif last_id_tmp == ():
            raise ValueError(f"Failed to retrieve value for last id: {retrieving_query}")
        else:
            last_id = last_id_tmp[0][0]

    else:
        raise ValueError(f"Failed check of {table_name} configuration")

    return last_id

def execute_query(query: str, table_name: str, data_dict: Dict, table_conf, metadata_params) -> Tuple[Any, Any]:
    """
    This function execute the query in the target database, if the query is an insert query it will return the id of the row inserted.
    If the query is an update query it will return the id of the row updated. If the query is duplicated it will retrieve the id of the row.

    Args:
        query (str): query to be executed in the target database
        table_name (str): table name to be used for the insert/update operation
        data_dict (Dict): dictionary containing key:values to be insert/update in the DB, used to retrieve the id of the row when the query is duplicated

    Raises:
        ValueError: raise an error when the query is not executed successfully

    Returns:
        int: id of the row inserted or updated
    """

    # Connecting to db
    conn = pymysql.connect(**metadata_params)
    cur  = conn.cursor()
    # Getting id name
    cur.execute(f"SHOW KEYS FROM {table_name} WHERE Key_name = 'PRIMARY'")
    id_name = cur.fetchone()[4]

    try:
        # Execute query and get id
        cur.execute(query)
        id_value = cur.lastrowid
    except pymysql.IntegrityError as e:
        if e.args[0] == 1062:
            # Retrieve value of a already exiting row
            logging.info(f"Query was duplicated in table {table_name}, retrieving id of the row")
            id_value = retrieve_row_id(data_dict, table_name, table_conf, metadata_params)

    except Exception as ee:
        print(f'Error: {ee}')
        raise ValueError from ee

    cur.close()
    conn.close()

    return id_value, id_name

def main():
    """Module's entry point"""
    logging.basicConfig(filename="write2db.log", level=logging.DEBUG, filemode='w',
                    format="%(asctime)s:%(levelname)s:%(message)s")
    parser = argparse.ArgumentParser(prog="write2db.py", 
                                    description="Create an insert or update queries and execute them in the target DB")
    parser.add_argument("--file-path", type=str,
                        help="Path to the JSON file containing data to insert or update in a DB")
    parser.add_argument('--update', action='store_true',
                        help="If this option is added it indicates the input data will be used to update")
    parser.add_argument('--empty', action='store_true',
                        help="If this option is added it indicates empty input data is allowed")
    parser.add_argument('--config', type=str,
                        help="Path to the JSON file containing the configuration of the tables")
    parser.add_argument('--metadata', type=str,
                        help="Path to the JSON file containing the metadata parameters")
    # Parsing arguments
    args = parser.parse_args()
    logging.info(f"Arguments: {args}")

    # Loading files
    logging.info(f"Loading file: {args.file_path}")
    file = open(str(args.file_path))
    input_data = json.load(file)
    file.close()

    if not input_data and args.empty:
        logging.info("Input data is empty. There is not data to process")
    elif not input_data and not args.empty:
        logging.info("Input data is empty. Data was expected in the file")
        raise ValueError("Input file is empty. Data expected or empty option was not provided")

    if args.config:
        if not os.path.exists(args.config):
            logging.info(f"File {args.config} does not exist")
            raise ValueError(f"File {args.config} does not exist")
        else:
            with open(args.config) as file:
                table_conf = json.load(file)

    if args.metadata:
        if not os.path.exists(args.metadata):
            logging.info(f"File {args.metadata} does not exist")
            raise ValueError(f"File {args.metadata} does not exist")
        else:
            with open(args.metadata) as file:
                metadata_params = json.load(file)

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
        if check:
            logging.info("Lists of dictionaries detected, processing each dictionary")
            for row in input_data[table_name]:
                query = create_query(row, table_name, update, table_conf, metadata_params)
                id_value, id_name= execute_query(query, table_name, input_data[table_name], table_conf, metadata_params)
                logging.info(f"Data was inserted in {table_name}. Last value of {id_name} is {id_value}")
                # saving last id in dict
                last_id_dict.update({id_name:id_value})
        else:
            logging.info("Regular dictionary detected, processing key:value pair values")
            # Data is a dictionary (This part is not tested yet)
            query = create_query(input_data[table_name], table_name, update, table_conf, metadata_params)
            id_value, id_name= execute_query(query,table_name, input_data[table_name], table_conf, metadata_params)
            logging.info(f"Data was inserted in {table_name}. Last value of {id_name} is {id_value}")
            # saving last id in dict
            last_id_dict.update({id_name:id_value})

    # saving output
    with open(output, 'a') as file:
        json.dump(last_id_dict, file)
    file.close()

if __name__ == "__main__":
    main()