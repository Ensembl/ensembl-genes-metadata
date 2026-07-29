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

import argparse
import json
import logging
import os
from typing import Any, Dict, Tuple

import pymysql  # type:ignore


def escape_value(value) -> str:
    """Escape a value for safe embedding inside a single-quoted SQL string literal."""
    return pymysql.converters.escape_string(str(value))


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
                if all(
                    isinstance(item, dict) for item in value
                ):  # Check if all items in the list are dictionaries
                    # print(f"Key '{key}' is linked to a list of dictionaries.")
                    dict_islist = True
                else:
                    raise ValueError(f"Table '{key}' is not properly formatted")
            else:
                # print(f"Key '{key}' is not linked to a list of dictionaries.")
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

    if update:  # check update key -> ukey
        key = "ukey"
    else:  # check dependent key -> dkey
        key = "dkey"

    # If insert query do not require dkey:
    if table_conf[table_name][key] == "None":
        logging.info("The %s table does not require any dependent/update key. Update %s ", table_name, update)

    # If insert query do require dkey:
    elif table_conf[table_name][key] != "None":
        if table_conf[table_name][key] in data_dict.keys():
            logging.info("Key %s exist in %s. Update %s", table_conf[table_name][key], table_name, update)
        else:
            raise ValueError(
                f"key {table_conf[table_name][key]} not found in {table_name} table. Update {update}"
            )
    else:
        raise ValueError(
            f"Unexpected value {table_conf[table_name][key]} for table {table_name} table. Update {update}"
        )


def insert_query(data_dict: Dict, table_name: str, table_conf, metadata_params) -> str:
    """This functions create a mysql insert queries using the input data provided

    Args:
        data_dict (dic): dictionary containing the key:values to be inserted
        table_name (str): table that is used to insert new data

    Returns:
        str: returns mysql query
    """

    if table_conf[table_name]["method"] not in ["per_row", "per_row_key"]:
        # crete basic query
        logging.info("Creating basic query for table %s", table_name)
        table_var_string = ", ".join(list(data_dict.keys()))
        values_strings = ",".join([f"'{escape_value(value)}'" for value in list(data_dict.values())]).replace(
            "''", "NULL"
        )
        return f"""INSERT INTO {table_name} ({table_var_string}) VALUES ({values_strings}) ;"""

    logging.info("%s is an attribute table (key:value pairs) ", table_name)

    dkey = table_conf[table_name]["dkey"]
    logging.info("The dkey is %s", dkey)

    # Getting columns names
    conn = pymysql.connect(**metadata_params)
    cur = conn.cursor()
    cur.execute(f"SHOW COLUMNS FROM {table_name}")
    table_columns = cur.fetchall()
    columns = [column[0] for column in table_columns if column[3] != "PRI"]
    columns_string = ", ".join(columns)
    #
    value_list = []
    if dkey is None or dkey == "None":
        logging.info("%s is a per_row table without dkey %s", table_name, dkey)
        for key, value in data_dict.items():
            value_item = f"('{escape_value(key)}', '{escape_value(value)}')"
            value_list.append(value_item)
            values_string = ", ".join(value_list)

        return f"""INSERT INTO {table_name} ({columns_string}) VALUES {values_string}"""

    logging.info("%s is a per_row table with a dkey %s", table_name, dkey)
    dkey_value = data_dict[dkey]
    for key, value in data_dict.items():
        if key != dkey:
            if table_conf[table_name]["method"] in ["per_row"]:
                value_item = f"('{escape_value(dkey_value)}', '{escape_value(key)}', '{escape_value(value)}')"
            elif table_conf[table_name]["method"] in ["per_row_key"]:
                logging.info("%s is an attribute table (key only) ", table_name)
                value_item = f"('{escape_value(dkey_value)}', '{escape_value(key)}')"
            else:
                raise ValueError(
                    f"Invalid value in table config - method: {table_conf[table_name]['method'] } "
                )

            value_list.append(value_item)
            values_string = ", ".join(value_list)
    return f"""INSERT INTO {table_name} ({columns_string}) VALUES {values_string}"""


def update_query(data_dict: Dict, table_name: str, table_conf) -> str:
    """
    This functions create a mysql update query using the input data provided

    Args:
        data_dict (dict): dictionary containing the key:values to be inserted
        table_name (str): table that is being updated

    Raises:
        ValueError: when the update key (ukey) is not present in data_dict

    Returns:
        str: returns a mysql query
    """
    update_list = []
    condition = None
    for key, value in data_dict.items():
        if table_conf[table_name]["ukey"] == key:
            condition = f"{key} = {value}"
        else:
            update_list.append(f"{key} = '{escape_value(value)}'")

    if condition is None:
        raise ValueError(f"Update key not found in the provided data for table {table_name}")

    update_values = ",".join(update_list)

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

    if update:  # input data will be used to update a row
        query = update_query(data_dict, table_name, table_conf)
    else:  # input data will be used to insert a new row
        query = insert_query(data_dict, table_name, table_conf, metadata_params)

    logging.info("Query created: %s", query)
    return query


def _retrieve_row_id_per_row_method(  # pylint: disable=too-many-branches
    data_dict, table_name, table_conf, cur, id_name
) -> int:
    """Retrieve (or insert-then-retrieve) the row id for a per_row / per_row_key table."""
    dkey = table_conf[table_name]["dkey"]

    cur.execute(f"SHOW COLUMNS FROM {table_name}")
    table_columns = cur.fetchall()
    columns = [column[0] for column in table_columns if column[3] != "PRI"]
    dkey_value = "None" if (dkey is None or dkey == "None") else data_dict[dkey]

    logging.info("Retriving IDs of inserted data. Table:%s. dkey: %s ", table_name, dkey_value)
    last_id = None
    for key, value in data_dict.items():
        if key == dkey:
            continue

        if table_conf[table_name]["method"] in ["per_row"] and dkey_value != "None":
            condition_string = (
                f"{columns[0]} = '{escape_value(dkey_value)}' AND {columns[1]} =  '{escape_value(key)}' "
                f"AND {columns[2]} = '{escape_value(value)}'"
            )
        elif table_conf[table_name]["method"] in ["per_row_key"]:
            condition_string = (
                f"{columns[0]} = '{escape_value(dkey_value)}' AND {columns[1]} =  '{escape_value(key)}'"
            )
        elif table_conf[table_name]["method"] in ["per_row"] and dkey_value == "None":
            condition_string = (
                f"{columns[0]} =  '{escape_value(key)}' AND {columns[1]} = '{escape_value(value)}'"
            )
        else:
            raise ValueError(f"Invalid value in table config - method: {table_conf[table_name]['method'] } ")

        retrieving_query = f"SELECT {id_name} FROM {table_name} WHERE {condition_string} ;"
        logging.info("Retriving IDs of inserted data. Query: %s", retrieving_query)
        cur.execute(retrieving_query)
        last_id_tmp = cur.fetchall()

        if len(last_id_tmp) > 1:
            raise ValueError(
                f"The query retrieves more than more value, unique value expected {retrieving_query}"
            )
        if last_id_tmp == () and table_conf[table_name]["method"] in ["per_row"]:
            logging.info("Failed to retrieve value for last id. Inserting missing data")
            if table_name == "taxonomy":
                query_missing_insert = (
                    f"UPDATE {table_name} SET {columns[1]} = '{escape_value(key)}' "
                    f"WHERE {columns[0]} = '{escape_value(dkey_value)}' "
                    f"AND {columns[2]} = '{escape_value(value)}' ;"
                )
            elif table_name == "taxonomy_name":
                query_missing_insert = (
                    f"INSERT INTO {table_name} ({columns[0]}, {columns[1]}) "
                    f"VALUES ('{escape_value(key)}', '{escape_value(value)}') ;"
                )
            else:
                query_missing_insert = (
                    f"INSERT INTO {table_name} ({columns[0]}, {columns[1]}, {columns[2]}) "
                    f"VALUES ('{escape_value(dkey_value)}', '{escape_value(key)}', '{escape_value(value)}') ;"
                )
            logging.info("Insert/update query: %s", query_missing_insert)
            cur.execute(query_missing_insert)
            last_id = cur.lastrowid
        else:
            logging.info("Retrieved value for last id: %s", last_id_tmp)
            last_id = last_id_tmp[0][0]

    if last_id is None:
        raise ValueError(f"Could not determine a last id for table {table_name}: no non-dkey values provided")

    return last_id


def _retrieve_row_id_per_col_method(data_dict, table_name, constraint, cur, id_name) -> int:
    """Retrieve the row id for a per_col table, using its uniqueness constraint."""
    condition_list = [f"{key[0]} = '{escape_value(data_dict[key[0]])}'" for key in constraint]
    condition_string = " AND ".join(condition_list)

    retrieving_query = f"SELECT {id_name} FROM {table_name} WHERE {condition_string} ;"
    cur.execute(retrieving_query)
    last_id_tmp = cur.fetchall()

    if len(last_id_tmp) > 1:
        raise ValueError(
            f"The query retrieves more than more value, unique value expected {retrieving_query}"
        )
    if last_id_tmp == ():
        raise ValueError(f"Failed to retrieve value for last id: {retrieving_query}")

    return last_id_tmp[0][0]


def retrieve_row_id(data_dict: Dict, table_name: str, table_conf, metadata_params) -> int:
    """
    If data is already inserted in the db, this function retrieves the id of the row using the
    data provided and the table name, it will retrieve the uniqueness constrain of the table and
    use it to retrieve the id of the row.

    Args:
        data_dict (dict): dictionary containing key:values to be insert/update in the DB
        table_name (str): table name to be used for the insert/update operation

    Returns:
        int: id of the row
    """
    logging.info("Retriving IDs of inserted data")
    # Establishing connection to DB
    conn = pymysql.connect(**metadata_params)
    cur = conn.cursor()

    # Getting id key name
    query = f"SHOW KEYS FROM {table_name} WHERE Key_name = 'PRIMARY'"
    cur.execute(query)
    id_name = cur.fetchone()[4]
    logging.info("Retrieving value of %s from %s", id_name, table_name)

    # Getting constraint keys
    constraint_query = f"""SELECT column_name FROM information_schema.key_column_usage
    WHERE
        table_schema = 'gb_assembly_metadata'
        AND table_name = '{table_name}'
        AND constraint_name != 'PRIMARY'
        AND referenced_table_name IS NULL;"""
    cur.execute(constraint_query)
    constraint = cur.fetchall()
    logging.info(" Detected uniqueness constrains: %s", constraint)

    method = table_conf[table_name]["method"]
    if method in ["per_row", "per_row_key"]:
        last_id = _retrieve_row_id_per_row_method(data_dict, table_name, table_conf, cur, id_name)
    elif method == "per_col":
        last_id = _retrieve_row_id_per_col_method(data_dict, table_name, constraint, cur, id_name)
    else:
        raise ValueError(f"Failed check of {table_name} configuration")

    cur.close()
    conn.close()

    return last_id


def execute_query(
    query: str, table_name: str, data_dict: Dict, table_conf, metadata_params
) -> Tuple[Any, Any]:
    """
    This function execute the query in the target database, if the query is an insert query it will
    return the id of the row inserted. If the query is an update query it will return the id of the
    row updated. If the query is duplicated it will retrieve the id of the row.

    Args:
        query (str): query to be executed in the target database
        table_name (str): table name to be used for the insert/update operation
        data_dict (Dict): dictionary containing key:values to be insert/update in the DB, used to
            retrieve the id of the row when the query is duplicated

    Raises:
        ValueError: raise an error when the query is not executed successfully

    Returns:
        int: id of the row inserted or updated
    """

    # Connecting to db
    conn = pymysql.connect(**metadata_params)
    cur = conn.cursor()
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
            logging.info("Query was duplicated in table %s, retrieving id of the row", table_name)
            id_value = retrieve_row_id(data_dict, table_name, table_conf, metadata_params)

    except Exception as ee:
        print(f"Error: {ee}")
        raise ValueError from ee

    cur.close()
    conn.close()

    return id_value, id_name


def load_input_data(args: argparse.Namespace) -> Tuple[Dict, Dict, Dict]:
    """Load and validate the input JSON, table config, and metadata DB parameters from CLI args."""
    with open(str(args.file_path), encoding="utf-8") as file:
        input_data = json.load(file)

    # Checking input file
    if not input_data and args.empty:
        logging.info("Input data is empty. There is not data to process")
    elif not input_data and not args.empty:
        logging.info("Input data is empty. Data was expected in the file")
        raise ValueError("Input file is empty. Data expected or empty option was not provided")

    # Check Configuration File
    if not args.config or not os.path.exists(args.config):
        raise ValueError(f"Please provide a valid --config file path, got: {args.config!r}")
    with open(args.config, encoding="utf-8") as file:
        table_conf = json.load(file)

    # Check metadata DB params
    if not args.metadata or not os.path.exists(args.metadata):
        raise ValueError(f"Please provide a valid --metadata file path, got: {args.metadata!r}")
    with open(args.metadata, encoding="utf-8") as file:
        metadata_params = json.load(file)

    return input_data, table_conf, metadata_params


def process_tables(input_data: Dict, update: bool, table_conf: Dict, metadata_params: Dict) -> Dict:
    """Create and execute an insert/update query for every table in input_data."""
    last_id_dict = {}
    # Process each key from input json file
    for table_name in input_data:
        logging.info("Processing input data, loading %s table", table_name)
        # Check input data structure
        check = check_dict_structure(input_data[table_name])
        if check:
            logging.info("Lists of dictionaries detected, processing each dictionary")
            for row in input_data[table_name]:
                query = create_query(row, table_name, update, table_conf, metadata_params)
                id_value, id_name = execute_query(
                    query, table_name, input_data[table_name], table_conf, metadata_params
                )
                logging.info("Data was inserted in %s. Last value of %s is %s", table_name, id_name, id_value)
                # saving last id in dict
                last_id_dict.update({id_name: id_value})
        else:
            logging.info("Regular dictionary detected, processing key:value pair values")
            # Data is a dictionary (This part is not tested yet)
            query = create_query(input_data[table_name], table_name, update, table_conf, metadata_params)
            id_value, id_name = execute_query(
                query, table_name, input_data[table_name], table_conf, metadata_params
            )
            logging.info("Data was inserted in %s. Last value of %s is %s", table_name, id_name, id_value)
            # saving last id in dict
            last_id_dict.update({id_name: id_value})

    return last_id_dict


def main():
    """Module's entry point"""
    logging.basicConfig(
        filename="write2db.log",
        level=logging.DEBUG,
        filemode="w",
        format="%(asctime)s:%(levelname)s:%(message)s",
    )
    parser = argparse.ArgumentParser(
        prog="write2db.py", description="Create an insert or update queries and execute them in the target DB"
    )
    parser.add_argument(
        "--file-path", type=str, help="Path to the JSON file containing data to insert or update in a DB"
    )
    parser.add_argument(
        "--update",
        action="store_true",
        help="If this option is added it indicates the input data will be used to update",
    )
    parser.add_argument(
        "--empty",
        action="store_true",
        help="If this option is added it indicates empty input data is allowed",
    )
    parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to the JSON file containing the configuration of the tables",
    )
    parser.add_argument(
        "--metadata",
        type=str,
        required=True,
        help="Path to the JSON file containing the metadata parameters",
    )
    # Parsing arguments
    args = parser.parse_args()
    logging.info("Arguments: %s", args)

    # Loading files
    logging.info("Loading file: %s", args.file_path)
    input_data, table_conf, metadata_params = load_input_data(args)

    # Get output name and create file
    root_name, _ = os.path.splitext(args.file_path)
    output = root_name + ".last_id"
    with open(output, "w", encoding="utf-8"):
        pass

    last_id_dict = process_tables(input_data, args.update, table_conf, metadata_params)

    # saving output
    with open(output, "a", encoding="utf-8") as file:
        json.dump(last_id_dict, file)


if __name__ == "__main__":
    main()
