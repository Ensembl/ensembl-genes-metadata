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

import logging
import argparse
import os
import pymysql  # type: ignore
from datetime import datetime
import json


def execute_query(query, db_params):
    conn = pymysql.connect(**db_params)
    cursor = conn.cursor()
    cursor.execute(query)
    result = cursor.fetchall()
    cursor.close()
    conn.close()
    return result


def fetch_gca_list(metadata_params, full_screen, screen_date):
    """Get a list of GCAs to check their status and metadata
    Args:
        metadata_params (dict): dictionary containing the connecting params for the gb_assembly_metadata database
        full_screen (boolean): if true it will update all the available records in the DB. Otherwise, it will update vertebrate or assemblies from the biodiversity projects
        screen_date (str): date to filter assemblies released after the given date (format YYYY-MM-DD)
    Returns:
        list: GCAs list
    """

    if full_screen:
        logging.info("Running in full screen mode")
        query_get_gca_list = """
        SELECT CONCAT(a.GCA_chain, '.', a.GCA_version) AS GCA
        FROM assembly a
        WHERE a.is_current = 'current';
        """
    else:
        logging.info("Running only high priority assemblies (vertebrates or relevant bioprojects)")
        query_get_gca_list = f"""
        SELECT DISTINCT CONCAT(a.GCA_chain, '.', a.GCA_version) AS GCA
        FROM bioproject b
        JOIN assembly a ON a.assembly_id = b.assembly_id
        JOIN species s ON a.lowest_taxon_id = s.lowest_taxon_id
        JOIN taxonomy t ON t.lowest_taxon_id = a.lowest_taxon_id
        WHERE s.scientific_name NOT LIKE '% sp.%'
            AND a.is_current = 'current'
            AND a.release_date > '{screen_date}'
            AND (t.taxon_class_id = '7711'
                OR b.bioproject_id IN (
                    'PRJNA533106', 'PRJEB40665', 'PRJEB61747',
                    'PRJEB43510', 'PRJEB47820', 'PRJNA813333', 'PRJNA489243'

        ));
        """

    # Get list from assembly_metadata and parse
    fetch_gca = execute_query(query_get_gca_list, metadata_params)
    gca_list = [gca[0] for gca in fetch_gca]
    logging.info(f"Total number of assemblies to check: {len(gca_list)}")

    return gca_list


def main():
    """module's entry-point"""

    logging.basicConfig(
        filename="fetch_assemblies.log",
        level=logging.DEBUG,
        filemode="w",
        format="%(asctime)s:%(levelname)s:%(message)s",
    )

    parser = argparse.ArgumentParser(
        prog="fetch_assemblies.py", description="Fetch a list of GCAs to run an update of metadata."
    )

    parser.add_argument("--metadata", type=str, help="Path to the metadata database params in json format")
    parser.add_argument(
        "--full_screen",
        action="store_true",
        help="If set, it will update all the available records in the DB",
    )
    parser.add_argument(
        "--screen_date",
        type=str,
        default="2019-01-01",
        help="If set, it will update all the records released after the given date (format YYYY-MM-DD). Default is 2019-01-01",
    )

    args = parser.parse_args()
    logging.info(args)

    if args.metadata:
        if not os.path.exists(args.metadata):
            raise ValueError("Please enter a valid file path for metadata database parameters")
        else:
            with open(args.metadata, "r") as file:
                metadata_params = json.load(file)

    if args.screen_date:
        if not datetime.strptime(args.screen_date, "%Y-%m-%d"):
            raise ValueError("Please enter a valid date format (YYYY-MM-DD)")
        else:
            logging.info(f"Custom date provided to retrieve assemblies: {args.screen_date}")
    else:
        logging.info("Default date will be used to retrieve assemblies: 2019-01-01")

    if args.full_screen:
        logging.info("Full screen mode activated")
    else:
        logging.info("Only high priority assemblies (vertebrates or relevant bioprojects) will be processed")

    gca_list = fetch_gca_list(metadata_params, args.full_screen, args.screen_date)

    if len(gca_list) > 0:

        with open(f"assemblies_to_update_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt", "w") as file:
            for gca in gca_list:
                print(gca)
                file.write(gca + "\n")
        file.close()

        logging.info(
            f"Accessions to register: {len(gca_list)}. Please note that some assemblies might belong to unspecified species"
        )

    else:
        logging.info(f"No assemblies found since {args.screen_date}")


if __name__ == "__main__":
    main()
