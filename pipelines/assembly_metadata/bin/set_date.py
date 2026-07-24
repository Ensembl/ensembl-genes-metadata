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

"""Determine the date to use when fetching new assemblies, based on the
update_date table in the assembly metadata database.

Args:
    metadata (str): path to the metadata database connection parameters in json format
    full_screen (bool): if set, retrieve the full-screen update date instead of the regular update date

Returns:
    stdout: the date to use, formatted MM/DD/YYYY
"""

import argparse
import json
import logging

import pymysql  # type: ignore


def get_full_screen_date(metadata_params: dict) -> str:
    """Retrieve the date recorded for the last full-screen update."""
    query = "SELECT DATE_FORMAT(date_value, '%m/%d/%Y') FROM update_date WHERE update_type = 'full_screen'"
    with pymysql.connect(**metadata_params) as conn:
        with conn.cursor() as cur:
            cur.execute(query)
            row = cur.fetchone()
    return row[0]


def get_regular_update_date(metadata_params: dict) -> str:
    """Retrieve the date to use for a regular update: one day before the last recorded update."""
    query = (
        "SELECT DATE_FORMAT(DATE_SUB(date_value, INTERVAL 1 DAY), '%m/%d/%Y') FROM update_date "
        "WHERE update_type = 'regular_update'"
    )
    with pymysql.connect(**metadata_params) as conn:
        with conn.cursor() as cur:
            cur.execute(query)
            row = cur.fetchone()
    return row[0]


def main():
    """Module's entry point"""
    logging.basicConfig(
        filename="set_date.log",
        level=logging.DEBUG,
        filemode="w",
        format="%(asctime)s:%(levelname)s:%(message)s",
    )

    parser = argparse.ArgumentParser(
        prog="set_date.py",
        description="Determine the date to use when fetching new assemblies.",
    )
    parser.add_argument(
        "--metadata", type=str, required=True, help="Path to the metadata database params in json format"
    )
    parser.add_argument(
        "--full_screen",
        action="store_true",
        help="If set, retrieve the full-screen update date instead of the regular update date",
    )

    args = parser.parse_args()
    logging.info(args)

    with open(args.metadata, "r", encoding="utf-8") as file:
        metadata_params = json.load(file)

    if args.full_screen:
        date_value = get_full_screen_date(metadata_params)
    else:
        date_value = get_regular_update_date(metadata_params)

    print(date_value)


if __name__ == "__main__":
    main()
