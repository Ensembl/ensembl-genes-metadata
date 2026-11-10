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


def comparing_refseq(data, accession, metadata_params):
    """Compare available paired REFSEQ from NCBI and Registry and generate update queries if needed.
    Args:
        data (dict): JSON data retrieved from NCBI API for the assembly
        accession (str): GCA accession of the assembly
        metadata_params (dict): Database connection parameters for metadata database
    Returns:
        str: A line with accession, current refseq in Registry, and refseq from NCBI
    """

    # Getting the paired RefSeq assembly and its status from NCBI.  The
    # top-level report is for the requested GCA assembly; paired_assembly is
    # the corresponding GCF record.
    report = data["reports"][0]
    assembly_info = report.get("assembly_info", {})
    paired_assembly = assembly_info.get("paired_assembly", {})
    paired_accession = report.get("paired_accession", "") or paired_assembly.get("accession", "")
    paired_status = paired_assembly.get("status", "unknown")

    # An assembly without a paired RefSeq accession in NCBI does not require
    # any change in Registry.
    if not paired_accession:
        logging.info(f"No update needed for RefSeq accession of assembly {accession}")
        return f"{accession}, refseq_check, false, NA, NA"

    # Getting the assembly ID and the most recently registered RefSeq accession
    # from Registry.  The latest row is the current row, regardless of whether
    # its GCF status is current, suppressed, or another NCBI status.
    # The RefSeq accession is now maintained in refseq_assembly rather than in
    # assembly.  A LEFT JOIN is needed for assemblies which have no RefSeq row.
    query_refseq = (
        "SELECT a.assembly_id, ra.refseq_accession, ra.is_current "
        "FROM assembly a "
        "LEFT JOIN refseq_assembly ra "
        "ON ra.assembly_id = a.assembly_id "
        f"WHERE CONCAT(a.gca_chain, '.', a.gca_version) = '{accession}' "
        "ORDER BY ra.refseq_assembly_id DESC LIMIT 1"
    )
    assembly_rows = execute_query(query_refseq, metadata_params)
    if not assembly_rows:
        raise ValueError(f"Assembly {accession} is not registered")

    assembly_id, existing_pair, existing_status = assembly_rows[0]

    # Comparison
    if existing_pair == paired_accession and existing_status == paired_status:
        logging.info(f"No update needed for RefSeq accession of assembly {accession}")
        output_line = f"{accession}, refseq_check, false, NA, NA"
    elif not existing_pair:
        logging.info(
            f"No RefSeq accession found for assembly {accession} in Registry, setting to {paired_accession} ({paired_status})"
        )
        query_insert_refseq = (
            "INSERT INTO refseq_assembly "
            f"(assembly_id, refseq_accession, is_current) VALUES ({assembly_id}, '{paired_accession}', '{paired_status}');"
        )
        logging.info(query_insert_refseq)
        execute_write(query_insert_refseq, metadata_params)
        output_line = f"{accession}, refseq_check, true, no_refseq, {paired_accession}"
    elif existing_pair == paired_accession:
        logging.info(
            f"Update needed for RefSeq status of assembly {accession}: current status is {existing_status}, GCF status is {paired_status}"
        )
        query_update_status = (
            "UPDATE refseq_assembly SET is_current = "
            f"'{paired_status}' WHERE assembly_id = {assembly_id} "
            f"AND refseq_accession = '{existing_pair}';"
        )
        logging.info(query_update_status)
        execute_write(query_update_status, metadata_params)
        output_line = f"{accession}, refseq_check, true, {existing_status}, {paired_status}"
    else:
        logging.info(
            f"Update needed for RefSeq accession of assembly {accession}: current RefSeq in Registry is {existing_pair}, RefSeq from NCBI is {paired_accession} ({paired_status})"
        )
        query_update_old = (
            "UPDATE refseq_assembly SET is_current = 'previous' "
            f"WHERE assembly_id = {assembly_id} AND refseq_accession = '{existing_pair}';"
        )
        query_insert_refseq = (
            "INSERT INTO refseq_assembly "
            f"(assembly_id, refseq_accession, is_current) VALUES ({assembly_id}, '{paired_accession}', '{paired_status}');"
        )
        logging.info(query_update_old)
        execute_write(query_update_old, metadata_params)
        logging.info(query_insert_refseq)
        execute_write(query_insert_refseq, metadata_params)
        output_line = f"{accession}, refseq_check, true, {existing_pair}, {paired_accession}"

    return output_line


def main():
    """Module's entry point"""

    logging.basicConfig(
        filename="update_assembly_refseq.log",
        level=logging.DEBUG,
        filemode="w",
        format="%(asctime)s:%(levelname)s:%(message)s",
    )

    parser = argparse.ArgumentParser(
        prog="update_assembly_refseq.py",
        description="Checks if there is an new available paired refseq accession.",
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

    args = parser.parse_args()
    logging.info(args)

    accession = args.accession.strip()

    with open(args.accession_json, "r") as json_file:
        data = json.load(json_file)

    metadata_params = args.metadata_params

    output_line = comparing_refseq(data, accession, metadata_params)
    print(output_line)


if __name__ == "__main__":
    main()
