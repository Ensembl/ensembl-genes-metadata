#!/usr/bin/env python3
# pylint: disable=missing-module-docstring
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
from typing import Dict, List, Union


def check_dict_structure(input_dict: Union[Dict, List[Dict]]) -> bool:
    """This function checks the structure of the database,
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
                if all(isinstance(item, dict) for item in value):
                    dict_list = True
                else:
                    raise ValueError(f"Table '{key}' is not properly formatted")
            else:
                dict_list = False
    return dict_list


def check_key(data_dict: Dict, table_name: str, update: bool, table_conf: Dict[str, Dict[str, str]]) -> None:
    """
    This function checks if the data has all keys necessary for a successful execution.
    It is used for insert queries or update queries

    Args:
        data_dict (dict): input data of one table to load in db
        table_name (str): table to load in db
        update (boolean): True if the data is to update an existing row. False is default (insert query)
        table_conf (dict): configuration dictionary containing keys for each table

    Raises:
        ValueError: raise an error when keys are missing
    """
    if update:
        key = "ukey"
    else:
        key = "dkey"

    if table_conf[table_name][key] == "None":
        logging.info("The %s table does not require any dependent/update key. Update %s ", table_name, update)
    elif table_conf[table_name][key] != "None":
        if table_conf[table_name][key] in data_dict.keys():
            logging.info("Key %s exists in %s. Update %s", table_conf[table_name][key], table_name, update)
        else:
            raise ValueError(
                f"key {table_conf[table_name][key]} not found in {table_name} table. Update {update}"
            )
    else:
        raise ValueError(
            f"Unexpected value {table_conf[table_name][key]} for table {table_name} table. Update {update}"
        )


def insert_query(data_dict: Dict, table_name: str) -> str:
    """This function creates a MySQL insert query using the input data provided

    Args:
        data_dict (dict): dictionary containing the key:values to be inserted
        table_name (_type_): table that is used to insert new data

    Returns:
        str: returns MySQL query
    """
    table_var_string = ", ".join(list(data_dict.keys()))
    values_strings = ",".join([f"'{value}'" for value in list(data_dict.values())]).replace("''", "NULL")

    return f"""INSERT IGNORE INTO {table_name} ({table_var_string}) VALUES ({values_strings}) ;"""


def update_query(data_dict: Dict, table_name: str, table_conf: Dict[str, Dict[str, str]]) -> str:
    """
    This function creates a MySQL update query using the input data provided

    Args:
        data_dict (dict): dictionary containing the key:values to be inserted
        table_name (str): table that is being updated
        table_conf (dict): configuration dictionary containing keys for each table

    Returns:
        str: returns a MySQL query
    """
    update_list = []
    for key, value in data_dict.items():
        if table_conf[table_name]["ukey"] != key:
            update_list.append(f"{key} = {value}")
        else:
            condition = f"{key} = {value}"

    update_values = ",".join(update_list)

    return f"UPDATE {table_name} SET {update_values} WHERE {condition}"


def create_query(
    data_dict: Dict, table_name: str, update: bool, table_conf: Dict[str, Dict[str, str]]
) -> str:
    """
    This function creates an insert or update MySQL query depending if update argument was provided (True).
    It uses the check_key function to determine if the data provided has enough keys to create the query.

    Args:
        data_dict (dict): dictionary containing key:values to be insert/update in the DB
        table_name (str): table name to be used for the insert/update operation
        update (bool): boolean variable indicating if the operation is an update
        table_conf (dict): configuration dictionary containing keys for each table

    Returns:
        str: MySQL query
    """
    check_key(data_dict, table_name, update, table_conf)

    if update:
        query = update_query(data_dict, table_name, table_conf)
    else:
        query = insert_query(data_dict, table_name)

    return query


def main():
    """Entry point"""

    table_conf = {
        "run": {"method": "per_col", "dkey": "None", "ukey": "None"},
        "study": {"method": "per_col", "dkey": "None", "ukey": "None"},
        "data_files": {"method": "per_col", "dkey": "run_id", "ukey": "file_id"},
        "align": {"method": "per_col", "dkey": "run_id", "ukey": "None"},
    }

    logging.basicConfig(
        filename="write2db.log",
        level=logging.DEBUG,
        filemode="w",
        format="%(asctime)s:%(levelname)s:%(message)s",
    )

    parser = argparse.ArgumentParser(prog="write.py", description="Write JSON files to DB")

    parser.add_argument(
        "--file_path", type=str, help="Path to the JSON file containing data to insert/update in a DB"
    )
    parser.add_argument(
                    "--output_dir", type=str, help="Output directory"
                        )
    parser.add_argument(
        "--update",
        type=str,
        choices=["True", "False"],
        default="False",
        help="Set to True if the input data will be used to update",
    )

    # Parsing arguments
    args = parser.parse_args()
    with open(str(args.file_path)) as input_file:
        input_data = json.load(input_file)

    update = args.update == "True"
    #with open(str(args.output_dir)+'/output_query.txt', 'w') as output_file:

        # Module
    for table_name in input_data:
            logging.info("Processing input data, loading %s table", table_name)
            # Check input data structure
            check = check_dict_structure(input_data[table_name])
            if check:
                logging.info("List of dictionaries detected, processing each dictionary")
                for row in input_data[table_name]:
                    query = create_query(row, table_name, update, table_conf)
                    print(query)
                    #output_file.write(query + '\n')  # Write query to file
            else:
                # Data is a dictionary
                query = create_query(input_data[table_name], table_name, update, table_conf)
                print(query)
                #output_file.write(query + '\n')  # Write query to file

    # Close the file after writing all queries
    #output_file.close()


if __name__ == "__main__":
    main()
