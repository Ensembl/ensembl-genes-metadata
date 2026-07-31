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


def taxonomy_checker(taxon_id: str, metadata_params: dict) -> bool:

    taxonomy_count = execute_query(
        f"SELECT COUNT(*) FROM taxonomy WHERE lowest_taxon_id = '{taxon_id}'",
        metadata_params,
    )[0][0]
    check = taxonomy_count == 7

    if check == False:
        query = f"DELETE FROM taxonomy WHERE lowest_taxon_id = '{taxon_id}'"
        execute_write(query, metadata_params)

    return check


def check_taxon_id(data, accession, metadata_params):

    # Taxon ID in NCBI
    taxon_id_ncbi = data["reports"][0].get("organism", "").get("tax_id")

    # Taxon ID in Registry
    query_taxon_id = (
        f"SELECT lowest_taxon_id from assembly WHERE CONCAT(gca_chain, '.', gca_version) = '{accession}'"
    )
    taxon_id = execute_query(query_taxon_id, metadata_params)[0][0]

    if taxon_id_ncbi == taxon_id:
        logging.info(f"No update needed for taxon ID of assembly {accession}")
        if taxonomy_checker(taxon_id, metadata_params):
            taxon_id_check = "pass"
        else:
            logging.info(f"Taxonomy hierarchy check failed for taxon ID of assembly {accession}")
            taxon_id_check = "fail"
    else:
        logging.info(
            f"Update needed for taxon ID of assembly {accession}: current taxon ID in Registry is {taxon_id}, taxon ID from NCBI is {taxon_id_ncbi}"
        )

        # Is the new taxon available in the species table
        query_new_taxon = f"SELECT COUNT(*) from species where lowest_taxon_id = {taxon_id_ncbi}"
        taxon_count = execute_query(query_new_taxon, metadata_params)[0][0]

        if taxon_count == 1:
            logging.info(f"New species taxon id exists in registry: {taxon_id_ncbi}")
            if taxonomy_checker(taxon_id_ncbi, metadata_params):
                taxon_id_check = "pass"
            else:
                taxon_id_check = "fail"

        elif taxon_count == 0:
            logging.info(
                f"New species taxon id not found in registry: {taxon_id_ncbi} --> Registry species and taxonomy hierarchy!"
            )
            taxon_id_check = "fail"
        else:
            raise ValueError(f"Species with taxon id {taxon_id_ncbi} detected multiple times: {taxon_count}")

    return ",".join([str(taxon_id), str(taxon_id_ncbi), taxon_id_check])


def comparing_basic_taxon_data(data, accession, metadata_params):

    output_line_list = []

    # Getting info from NCBI report
    taxon_id_ncbi = data["reports"][0].get("organism", "").get("tax_id")
    organism_name_ncbi = data["reports"][0].get("organism", "").get("organism_name")
    common_name_ncbi = data["reports"][0].get("organism", "").get("common_name", "")

    # Taxon ID check
    query_taxon_id = (
        f"SELECT lowest_taxon_id from assembly WHERE CONCAT(gca_chain, '.', gca_version) = '{accession}'"
    )
    taxon_id = execute_query(query_taxon_id, metadata_params)[0][0]

    if taxon_id_ncbi == taxon_id:
        logging.info(f"No update needed for taxon ID of assembly {accession}")
    else:
        logging.info(
            f"Update needed for taxon ID of assembly {accession}: current taxon ID in Registry is {taxon_id}, taxon ID from NCBI is {taxon_id_ncbi}"
        )
        query_update_taxon_id = f"""UPDATE assembly 
        SET lowest_taxon_id = '{taxon_id_ncbi}' 
        WHERE CONCAT(gca_chain, '.', gca_version) = '{accession}';"""
        logging.info(query_update_taxon_id)
        affected = execute_write(query_update_taxon_id, metadata_params)
        output_line = f"{accession}, taxon_id_check, true, {taxon_id}, {taxon_id_ncbi}"
        output_line_list.append(output_line)

    # Scientific and common name check
    query_species = (
        f"SELECT scientific_name, common_name from species where lowest_taxon_id = '{taxon_id_ncbi}';"
    )
    scientific_name, common_name = execute_query(query_species, metadata_params)[0]

    if scientific_name != organism_name_ncbi:
        logging.info(
            f"Update required for scientific name of taxon ID {taxon_id_ncbi}: current scientific name in Registry is {scientific_name}, scientific name from NCBI is {organism_name_ncbi}"
        )
        query_update_scientific_name = f"""UPDATE species 
        SET scientific_name = '{organism_name_ncbi}' 
        WHERE lowest_taxon_id = '{taxon_id_ncbi}';"""
        logging.info(query_update_scientific_name)
        affected = execute_write(query_update_scientific_name, metadata_params)
        output_line = f"{accession}, scientific_name_check, true, {scientific_name}, {organism_name_ncbi}"
        output_line_list.append(output_line)
    else:
        logging.info(f"No update required for scientific name of taxon ID {taxon_id_ncbi}")

    if common_name != common_name_ncbi:
        logging.info(
            f"Update required for common name of taxon ID {taxon_id_ncbi}: current common name in Registry is {common_name}, common name from NCBI is {common_name_ncbi}"
        )
        query_update_common_name = f"""UPDATE species 
        SET common_name = '{common_name_ncbi.replace("'", "''")}' 
        WHERE lowest_taxon_id = '{taxon_id_ncbi}';"""
        logging.info(query_update_common_name)
        affected = execute_write(query_update_common_name, metadata_params)
        output_line = f"{accession}, common_name_check, true, {common_name}, {common_name_ncbi}"
        output_line_list.append(output_line)
    else:
        logging.info(f"No update required for common name of taxon ID {taxon_id_ncbi}")

    return output_line_list


def main():
    """Module's entry point"""

    logging.basicConfig(
        filename="update_taxonomy.log",
        level=logging.DEBUG,
        filemode="w",
        format="%(asctime)s:%(levelname)s:%(message)s",
    )

    parser = argparse.ArgumentParser(
        prog="taxonomy.py",
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
        type=json.loads,
        required=True,
        help="JSON string with metadata database connection parameters",
    )
    parser.add_argument(
        "--taxonomy_check",
        action="store_true",
        help="When added the module will check if the taxon id of the assembly is recorded in the registry.",
    )
    parser.add_argument(
        "--taxonomy_update",
        action="store_true",
        help="When added the module will attempt to update the basic taxonomy information.",
    )

    args = parser.parse_args()
    logging.info(args)

    accession = args.accession.strip()

    with open(args.accession_json, "r") as json_file:
        data = json.load(json_file)

    metadata_params = args.metadata_params

    if not (args.taxonomy_check or args.taxonomy_update):
        raise ValueError("Select at least one mode: --taxonomy_check and/or --taxonomy_update.")

    # Taxonomy check
    if args.taxonomy_check:
        output_line = check_taxon_id(data, accession, metadata_params)
        print(output_line)

    if args.taxonomy_update:
        output_line_list = comparing_basic_taxon_data(data, accession, metadata_params)
        for output_line in output_line_list:
            print(output_line)


if __name__ == "__main__":
    main()
