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

import argparse
import json
import logging

from gb_metadata.db_utils import execute_query, execute_write


def comparing_status(data, accession, metadata_params):
    """Compare assembly status from NCBI and Registry and generate update queries if needed.
    Args:
        data (dict): JSON data retrieved from NCBI API for the assembly
        accession (str): GCA accession of the assembly
        metadata_params (dict): Database connection parameters for metadata database
    Returns:
        str: A line with accession, current status in Registry, and status from NCBI
    """

    # Getting info from Registry
    query_status = (
        f"SELECT is_current FROM assembly WHERE CONCAT(gca_chain, '.', gca_version) = '{accession}'"
    )
    is_current = execute_query(query_status, metadata_params)[0][0]

    # Getting info from NCBI
    assembly_status = data["reports"][0].get("assembly_info").get("assembly_status")
    warning = data["reports"][0].get("assembly_info").get("atypical", {}).get("warnings", "NA")

    # Comparing and generating update queries
    if assembly_status == is_current and warning == "NA":
        logging.info(
            f"No update needed for assembly {accession}. Current status in Registry and NCBI is {is_current}"
        )
        output_line = f"{accession}, asm_status, false, NA, NA"
    elif warning != "NA":
        logging.info(f"Assembly {accession} has warnings: {warning}. Update status based on warning.")
        query_update_status = f"UPDATE assembly a SET is_current = '{warning[0]}' WHERE CONCAT(a.gca_chain, '.', a.gca_version) = '{accession}';"
        logging.info(query_update_status)
        affected = execute_write(query_update_status, metadata_params)
        output_line = f"{accession}, asm_status, true, {is_current}, {warning}"
    else:
        logging.info(
            f"Update needed for assembly {accession}: current status in Registry is {is_current}, status from NCBI is {assembly_status}"
        )
        query_update_status = f"UPDATE assembly a SET is_current = '{assembly_status}' WHERE CONCAT(a.gca_chain, '.', a.gca_version) = '{accession}';"
        logging.info(query_update_status)
        affected = execute_write(query_update_status, metadata_params)
        output_line = f"{accession}, asm_status, true, {is_current}, {assembly_status}"

    return output_line


def main():
    """Module's entry point"""

    logging.basicConfig(
        filename="update_assembly_status.log",
        level=logging.DEBUG,
        filemode="w",
        format="%(asctime)s:%(levelname)s:%(message)s",
    )

    parser = argparse.ArgumentParser(
        prog="update_assembly_status.py",
        description="Retrieve metadata from NCBI API for a given GCA accession and store it in JSON files to be inserted in the database.",
    )

    parser.add_argument(
        "--accession_json",
        type=str,
        required=True,
        help="Path to JSON file with assembly metadata retrieved from NCBI API",
    )
    parser.add_argument("--accession", type=str, required=True, help="GCA accession to retrieve metadata")
    parser.add_argument(
        "--metadata_params",
        type=str,
        required=True,
        help="Database connection parameters for metadata database in JSON format",
    )

    args = parser.parse_args()
    logging.info(args)

    accession = args.accession.strip()

    with open(args.accession_json, "r") as json_file:
        data = json.load(json_file)

    with open(args.metadata_params, "r") as params_file:
        metadata_params = json.load(params_file)

    output_line = comparing_status(data, accession, metadata_params)
    print(output_line)


if __name__ == "__main__":
    main()
