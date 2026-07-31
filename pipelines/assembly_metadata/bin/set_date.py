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

from gb_metadata.db_utils import fetch_one_row


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
        "--metadata", type=json.loads, required=True, help="JSON string with metadata database connection parameters"
    )
    parser.add_argument(
        "--full_screen",
        action="store_true",
        help="If set, retrieve the full-screen update date",
    )

    args = parser.parse_args()
    logging.info(args)

    metadata_params = args.metadata

    if args.full_screen:
        query_full_screen = (
            "SELECT DATE_FORMAT(date_value, '%m/%d/%Y') " "FROM update_date WHERE update_type = 'full_screen'"
        )
        date_value = fetch_one_row(query_full_screen, metadata_params, "full_screen update date")[0]
    else:
        query_regular_date = (
            "SELECT DATE_FORMAT(DATE_SUB(date_value, INTERVAL 1 DAY), '%m/%d/%Y') "
            "FROM update_date WHERE update_type = 'regular_update'"
        )
        date_value = fetch_one_row(query_regular_date, metadata_params, "regular update date")[0]

    print(date_value)


if __name__ == "__main__":
    main()
