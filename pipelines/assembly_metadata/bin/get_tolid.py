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

from gb_metadata.db_utils import fetch_one_row
from gb_metadata.utils import connection_api


def get_tolid(taxon: int) -> str:
    """Retrieves the ToLID prefix for a given taxon from the ToLID API

    Args:
        taxon (int): lowest taxon id of the assembly

    Returns:
        str: ToLID prefix for the given taxon
    """

    logging.info(f"Connecting to ToLID API to retrieve prefix for taxon {taxon}")
    uri = f"https://id.tol.sanger.ac.uk/api/v2/species?taxonomyId={taxon}"
    response = connection_api(uri)
    data = response.json()
    logging.info(f"Response from ToLID API: {data}")

    if data["totalNumSpecies"] == 0:
        logging.info(f"No species found for taxon {taxon}")
        tolid_prefix = ""
    else:
        tolid_prefix = data["species"][0]["prefix"]

    return tolid_prefix


def main():
    """Module's entry point"""

    logging.basicConfig(
        filename="get_tolid.log",
        level=logging.DEBUG,
        filemode="w",
        format="%(asctime)s:%(levelname)s:%(message)s",
    )

    parser = argparse.ArgumentParser(
        prog="get_tolid.py", description="Get ToLID prefix for a given assembly based on the taxon ID"
    )
    parser.add_argument("--accession", type=str, help="Full GCA assembly accession")

    parser.add_argument(
        "--metadata",
        type=json.loads,
        required=True,
        help="JSON string with metadata database connection parameters",
    )

    args = parser.parse_args()
    logging.info(args)

    metadata_params = args.metadata

    logging.info(
        f"Connecting to metadata database to retrieve lowest_taxon_id and assembly_id for {args.accession}"
    )
    gca_chain = args.accession.split(".")[0]
    gca_version = args.accession.split(".")[1]
    query_taxon = f'SELECT assembly.assembly_id, species.species_taxon_id FROM assembly INNER JOIN species ON species.lowest_taxon_id = assembly.lowest_taxon_id WHERE gca_chain = "{gca_chain}" AND gca_version = "{gca_version}"'
    assembly_id, taxon = fetch_one_row(query_taxon, metadata_params, context=args.accession)

    tolid = get_tolid(taxon)

    logging.info(f"Connecting to metadata database to retrieve organism_id for {args.accession}")
    query_organism = f'SELECT organism_id FROM organism WHERE assembly_id = "{assembly_id}"'
    organism_id = fetch_one_row(query_organism, metadata_params, context=f"organism_id of {args.accession}")

    dtol_json = {"organism": {"organism_id": organism_id[0], "tol_prefix": tolid}}

    logging.info(f"Saving tolid json for {args.accession}")
    with open(f"{args.accession}_tolid.json", "w") as file:
        json.dump(dtol_json, file)
    file.close()


if __name__ == "__main__":
    main()
