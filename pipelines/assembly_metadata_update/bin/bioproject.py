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
from typing import Dict, Any
import pymysql  # type: ignore


def execute_query(query, db_params):
    conn = pymysql.connect(**db_params)
    cursor = conn.cursor()
    cursor.execute(query)
    result = cursor.fetchall()
    cursor.close()
    conn.close()
    return result


def execute_write(query: str, db_params: Dict[str, Any]) -> int:
    """
    Execute INSERT/UPDATE/DELETE and commit.
    Returns number of affected rows.
    """
    conn = pymysql.connect(**db_params)
    try:
        with conn.cursor() as cursor:
            affected = cursor.execute(query)
        conn.commit()
        return affected
    finally:
        conn.close()


def comparing_bioproject(data, accession, metadata_params):
    """Compare bioprojects from NCBI and Registry and generate insert queries if needed.
    Args:
        data (dict): JSON data retrieved from NCBI API for the assembly
        accession (str): GCA accession of the assembly
        metadata_params (dict): Database connection parameters for metadata database
    Returns:
        str: A line with updating information
    """

    # Getting info from Registry
    query_bioprojects = f"""SELECT bio.bioproject_id from assembly a 
    JOIN bioproject bio ON a.assembly_id = bio.assembly_id 
    WHERE CONCAT(a.gca_chain, '.', a.gca_version) = '{accession}'"""
    logging.info(query_bioprojects)
    bioprojects_registry = execute_query(query_bioprojects, metadata_params)
    logging.info(bioprojects_registry)

    bioproject_list = []
    for item in bioprojects_registry:
        bioproject_list.append(item[0])

    query_assembly_id = (
        f"SELECT assembly_id FROM assembly WHERE CONCAT(gca_chain, '.', gca_version) = '{accession}'"
    )
    logging.info(query_assembly_id)
    assembly_id = execute_query(query_assembly_id, metadata_params)[0][0]
    logging.info(assembly_id)

    # Getting info from NCBI
    bioproject_lineage = {}
    seen_accessions = set()

    bioproject_dict = data["reports"][0]["assembly_info"]["bioproject_lineage"][0]["bioprojects"]
    for item in bioproject_dict:
        bio_accession = item["accession"]
        title = item["title"].replace("'", "")
        # Check if accession is not seen before
        if bio_accession not in seen_accessions:
            seen_accessions.add(bio_accession)
            bioproject_lineage[bio_accession] = title

    # Comparing and generating update queries

    missing_bioprojects = list(set(bioproject_lineage) - set(bioproject_list))
    logging.info(missing_bioprojects)

    bioproject = []
    if missing_bioprojects:
        for bioproject in missing_bioprojects:
            title = bioproject_lineage[bioproject]
            logging.info(f"Insert bioproject {bioproject} - {title} for assembly {accession}")
            query_insert_bioproject = f"""INSERT INTO bioproject (assembly_id, bioproject_id) 
            VALUES ('{assembly_id}', '{bioproject}');"""
            logging.info(query_insert_bioproject)
            affected = execute_write(query_insert_bioproject, metadata_params)

        bioproject_line = "-".join(missing_bioprojects)
        output_line = f"{accession}, bioproject_check, false, NA, {bioproject_line}"
    else:
        logging.info(f"No bioproject insert needed for assembly {accession}")
        output_line = f"{accession}, bioproject_check, false, NA, NA"

    return output_line


def main():
    """Module's entry point"""

    logging.basicConfig(
        filename="update_bioproject.log",
        level=logging.DEBUG,
        filemode="w",
        format="%(asctime)s:%(levelname)s:%(message)s",
    )

    parser = argparse.ArgumentParser(
        prog="bioproject.py", description="Checks if the GCA has been added to other bioproject."
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

    output_line = comparing_bioproject(data, accession, metadata_params)
    print(output_line)


if __name__ == "__main__":
    main()
