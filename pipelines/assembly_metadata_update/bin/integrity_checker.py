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
import pymysql  # type:ignore
import logging
import json
import os


def execute_query(query, db_params):
    conn = pymysql.connect(**db_params)
    cursor = conn.cursor()
    cursor.execute(query)
    result = cursor.fetchall()
    cursor.close()
    conn.close()
    return result


def delete_assembly(assembly_id: int, metadata_params: dict):
    """
    Delete records from metadata database for a given assembly_id.
    Affects assembly, assembly_metrics, organism and bioproject tables

    Args:
        assembly_id (int): Assembly ID to be deleted
    """
    con = pymysql.connect(**metadata_params)
    cur = con.cursor()

    # Delete records from all metadata tables
    cur.execute(f"DELETE FROM assembly WHERE assembly_id = '{assembly_id}'")
    cur.execute(f"DELETE FROM assembly_metrics WHERE assembly_id = '{assembly_id}'")
    cur.execute(f"DELETE FROM organism WHERE assembly_id = '{assembly_id}'")
    cur.execute(f"DELETE FROM bioproject WHERE assembly_id = '{assembly_id}'")
    logging.info(f"Records for assembly_id {assembly_id} were deleted")
    con.close()


def records_checker(accession: str, metadata_params: dict, delete_records: bool = False) -> str:
    """
    Check whether an accession has the expected metadata records.

    Returns:
        True  -> all required records exist
        False -> one or more required records are missing (and may be deleted if delete_records=True)
    """
    logging.info(
        "Connecting to Assembly metadata database to retrieve the assembly_id of GCA: %s",
        accession,
    )

    query_assembly_id = (
        f"SELECT assembly_id FROM assembly WHERE CONCAT(gca_chain, '.', gca_version) = '{accession}';"
    )

    logging.info(query_assembly_id)
    rows = execute_query(query_assembly_id, metadata_params)

    if not rows:
        msg = f"Accession {accession} not found in assembly table. Nothing to check."
        logging.warning(msg)
        raise ValueError(msg)

    assembly_id = rows[0][0]

    query_lowest_taxon_id = (
        "SELECT lowest_taxon_id "
        "FROM assembly "
        f"WHERE CONCAT(gca_chain, '.', gca_version) = '{accession}'"
    )

    lowest_taxon_id = execute_query(query_lowest_taxon_id, metadata_params)[0][0]

    logging.info(
        "Checking records for assembly_id %s in assembly_metrics, organism and bioproject tables",
        assembly_id,
    )

    metrics_count = execute_query(
        f"SELECT COUNT(*) FROM assembly_metrics WHERE assembly_id = '{assembly_id}'",
        metadata_params,
    )[0][0]
    metrics_ok = metrics_count >= 1

    organism_count = execute_query(
        f"SELECT COUNT(*) FROM organism WHERE assembly_id = '{assembly_id}'",
        metadata_params,
    )[0][0]
    organism_ok = organism_count == 1

    bioproject_count = execute_query(
        f"SELECT COUNT(*) FROM bioproject WHERE assembly_id = '{assembly_id}'",
        metadata_params,
    )[0][0]
    bioproject_ok = bioproject_count >= 1

    logging.info(
        "Checking records for lowest_taxon_id %s in species and taxonomy table",
        assembly_id,
    )

    taxonomy_count = execute_query(
        f"SELECT COUNT(*) FROM taxonomy WHERE lowest_taxon_id = '{lowest_taxon_id}'",
        metadata_params,
    )[0][0]
    taxonomy_ok = taxonomy_count >= 1

    species_count = execute_query(
        f"SELECT COUNT(*) FROM species WHERE lowest_taxon_id = '{lowest_taxon_id}'",
        metadata_params,
    )[0][0]
    species_ok = species_count == 1

    all_ok = metrics_ok and organism_ok and bioproject_ok and taxonomy_ok and species_ok

    # Taxonomy names records
    taxonomy_name_count = execute_query(
        f"SELECT COUNT(*) FROM taxonomy_name WHERE taxon_class_id = '{lowest_taxon_id}'",
        metadata_params,
    )[0][0]

    if all_ok:
        status = "correct"

        if taxonomy_name_count == 0:
            status = "taxonomy_update"
            accession = lowest_taxon_id

    if not all_ok:
        missing = []
        if not metrics_ok:
            missing.append("assembly_metrics")
        if not organism_ok:
            missing.append("organism")
        if not bioproject_ok:
            missing.append("bioproject")
        if not taxonomy_ok:
            missing.append("taxonomy")
        if not species_ok:
            missing.append("species")

        genebuild_count = execute_query(
            f"SELECT COUNT(*) FROM genebuild_status where gca_accession = '{accession}'",
            metadata_params,
        )[0][0]

        if genebuild_count == 0:

            status = "delete"
            logging.info(f"Accession {accession} has missing data in {', '.join(missing)}.")

            if delete_records:
                delete_assembly(assembly_id, metadata_params)

        else:
            status = "check"
            logging.info(
                f"Accession {accession} has missing data in {', '.join(missing)}. But there is an annotation records, review manually."
            )

    return status, accession


def main():
    """
    Module's entry point
    """

    logging.basicConfig(
        filename="clean_gca_records.log",
        level=logging.DEBUG,
        filemode="w",
        format="%(asctime)s:%(levelname)s:%(message)s",
    )

    parser = argparse.ArgumentParser(
        prog="clean_gca_records.py",
        description="Clean records from metadata database based on GCA records completeness o a list of GCA accessions",
    )
    parser.add_argument("--accession", help="GCA accession")
    parser.add_argument(
        "--delete", default=False, action="store_true", help="Delete records with missing data"
    )
    parser.add_argument(
        "--metadata",
        help="Path to metadata database connection parameters",
    )

    args = parser.parse_args()
    logging.info(args)

    if args.metadata:
        if not os.path.exists(args.metadata):
            raise ValueError("Metadata params json file does not exist")
        else:
            with open(args.metadata, "r") as f:
                metadata_params = json.load(f)
                f.close()

    status, accession = records_checker(args.accession, metadata_params, args.delete)

    print(f"{status},{accession}")


if __name__ == "__main__":
    main()
