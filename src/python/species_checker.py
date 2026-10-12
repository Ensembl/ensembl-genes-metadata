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
The module add species_taxon_id and parlance_name to the species main key to json-like (.tmp) species file.
It retrieves the taxonomy classification and names to be store in the taxonomy and taxonomy_name tables.
The module when used in the pipeline requires --json-path, --ncbi_url and --enscode.
Additionally, to be used as standalone script to update the taxonomy tables, it requires
--taxonomy_update, --taxon_id, and --ncbi_url.

Raises:
    ValueError: invalid taxonomy rank
    ValueError: multiple prefix detected

Returns:
    str: a json file with the species information
"""

import argparse
import json
import logging
import os

from gb_metadata.utils import connection_api


def get_taxon_data(taxon_id: int, ncbi_url: str) -> dict:
    """It connects to the NCBI API to retrieve the taxonomy data of the lowest taxon id.

    Args:
        taxon_id (int): lowest taxon id of the assembly
        ncbi_url (str): valid url to connect to NCBI API, taxonomy endpoint

    Returns:
        dict: response from NCBI API
    """
    uri = f"{ncbi_url}/taxonomy/taxon/{taxon_id}/dataset_report"
    logging.info("URI: %s", uri)
    response = connection_api(uri)
    taxon_data = response.json()

    return taxon_data


def get_taxon_classification(taxon_data) -> tuple[dict, dict]:
    """
    This function retrieves the taxonomy classification from the NCBI API response.
    It returns the classification dictionary that will be user to update the taxonomy table.

    Args:
        taxon_data (dict): NCBI API response containing taxonomy information for the lowest taxon id.

    Returns:
        tuple[dict, dict]: classification dictionary and classification name dictionary
    """
    classification = taxon_data["reports"][0]["taxonomy"]["classification"]

    classification_dic = {}
    classification_name_dic = {}
    for rank in classification:
        if rank not in ["domain", "superkingdom"]:
            classification_dic.update({classification[rank]["id"]: rank})
            classification_name_dic.update(
                {classification[rank]["id"]: classification[rank]["name"].replace("'", "''")}
            )

    return classification_dic, classification_name_dic


def species_taxon(taxon_data, taxon_id) -> tuple[int, bool]:
    """
    This functions checks if the lowest taxon id is a species or a infraspecific taxon rank.
    It returns the species taxon id to be used in the species table.

    Args:
        lowest_taxon_id (str): the taxonomy id of the assembly, it is obtained from the assembly NCBI report

    Raises:
        ValueError: when the provided value is not a valid taxonomy rank. Valid values are
            species taxon ID or infraspecific taxon ID

    Returns:
        str: species taxon id, it could be the same lowest taxon id value
    """

    taxon_exists = True
    try:
        taxonomy = taxon_data["reports"][0]["taxonomy"]["rank"]
        if taxonomy in [
            "SUBSPECIES",
            "STRAIN",
            "VARIETAS",
            "GENOTYPE",
            "ISOLATE",
            "FORMA",
            "FORMA_SPECIALIS",
            "CLADE",
        ]:
            species_taxon_id = taxon_data["reports"][0]["taxonomy"]["classification"]["species"]["id"]
            logging.info("The assembly is a infraspecific taxon %s", taxonomy)
        elif taxonomy == "SPECIES":
            species_taxon_id = taxon_id
            logging.info("The assembly is a species taxon rank ")
        else:
            raise ValueError(f"Incorrect taxonomy ({taxonomy})")
    except KeyError as exc:
        if "errors" in taxon_data["reports"][0]:
            species_taxon_id = 0  # Set species taxon as zero to be identified by the reporting module
            taxon_exists = False
            logging.info("Taxon %s is not a recognized NCBI Taxonomy name", taxon_id)
        elif "taxonomy" in taxon_data["reports"][0]:
            logging.info(
                "Taxon do not have Rank available, retrieving information from another section of the report"
            )
            species_taxon_id = taxon_data["reports"][0]["taxonomy"]["classification"]["species"]["id"]
        else:
            raise KeyError(f"Taxon {taxon_id} retrieves an unexpected report") from exc

    return species_taxon_id, taxon_exists


def get_parlance_name(sci_name: str, enscode) -> str:
    """
    Search in the snp_static.txt file from the core_meta_update repository the parlance name
    for the scientific name provided

    Args:
        sci_name (str): scientific name of the species
        enscode (str): ENSCODE variable path

    Returns:
        str: parlance name when available
    """

    parlance_file = f"{enscode}/ensembl-genes/src/python/ensembl/genes/metadata/snp_static.txt"
    data_dict = {}

    logging.info("Reading parlance name file (snp_static.txt) to look for a match")
    with open(parlance_file, "r", encoding="utf-8") as file:
        for line in file:
            key, value = line.rsplit("\t", 1)
            data_dict[key.strip()] = value.strip()

    parlance_name = data_dict.get(sci_name, "")

    return parlance_name


def update_species_metadata(args) -> None:
    """Update the species JSON file with taxonomy classification and parlance name."""
    with open(args.json_path, "r", encoding="utf-8") as file:
        species_dict = json.load(file)

    if not args.enscode:
        raise ValueError("Please enter a valid path for ENSCODE")

    logging.info("Getting key values for the species: %s", species_dict["species"]["scientific_name"])
    # Get taxon data from NCBI API
    taxon_data = get_taxon_data(species_dict["species"]["lowest_taxon_id"], args.ncbi_url)
    species_taxon_id, taxon_exists = species_taxon(taxon_data, species_dict["species"]["lowest_taxon_id"])
    if taxon_exists:
        parlance_name = get_parlance_name(species_dict["species"]["scientific_name"], args.enscode)
        species_prefix = ""
        taxon_classification, _ = get_taxon_classification(taxon_data)
        taxon_classification_check = True
    else:
        logging.info(
            "Taxon do not exist in taxonomy: invalid lowest taxon id or assembly should be suppressed"
        )
        logging.info("Setting values to NA/NULL to later be detected by the integrity check")
        parlance_name = ""
        species_prefix = ""
        taxon_classification = {}
        taxon_classification_check = False

    # Update species dictionary with new values
    logging.info("Updating keys for species table")
    species_dict["species"].update(
        {
            "species_taxon_id": species_taxon_id,
            "parlance_name": parlance_name,
            "species_prefix": species_prefix,
        }
    )
    if taxon_classification_check:
        species_dict["taxonomy"] = taxon_classification
        species_dict["taxonomy"].update({"lowest_taxon_id": species_dict["species"]["lowest_taxon_id"]})

    # Saving results
    output_file = os.path.basename(args.json_path).replace(".tmp", ".json")
    logging.info("Saving output: %s", output_file)
    with open(output_file, "w", encoding="utf-8") as file:
        json.dump(species_dict, file)


def update_taxonomy_tables(args) -> None:
    """Update the taxonomy and taxonomy_name tables, used as a standalone script."""
    logging.info("Updating taxonomy table")
    taxon_dict = {}
    taxon_data = get_taxon_data(args.taxon_id, args.ncbi_url)
    taxon_classification, taxon_name_classification = get_taxon_classification(taxon_data)
    taxon_dict["taxonomy"] = taxon_classification
    taxon_dict["taxonomy"].update({"lowest_taxon_id": args.taxon_id})
    logging.info("Updating taxonomy name table")
    taxon_dict["taxonomy_name"] = taxon_name_classification

    species_taxon_id, _ = species_taxon(taxon_data, args.taxon_id)

    taxon_dict["species"] = {
        "lowest_taxon_id": args.taxon_id,
        "species_taxon_id": species_taxon_id,
        "scientific_name": taxon_data["reports"][0]["taxonomy"]["current_scientific_name"]["name"],
        "common_name": taxon_data["reports"][0]["taxonomy"].get("curator_common_name", ""),
    }

    # Saving results
    output_file_taxon = f"taxonomy_{args.taxon_id}.json"
    logging.info("Saving output: %s", output_file_taxon)
    with open(output_file_taxon, "w", encoding="utf-8") as file:
        json.dump(taxon_dict, file)


def main():
    """Module's entry point."""
    logging.basicConfig(
        filename="species_checker.log",
        level=logging.DEBUG,
        filemode="w",
        format="%(asctime)s:%(levelname)s:%(message)s",
    )
    parser = argparse.ArgumentParser(
        prog="species_checker.py", description="Update species related metadata."
    )
    parser.add_argument("--json-path", type=str, help="Path to the JSON-like (.tmp) species file")
    parser.add_argument("--ncbi_url", type=str, required=True, help="NCBI API URL")
    parser.add_argument("--enscode", type=str, help="ENSCODE path")
    parser.add_argument("--taxon_id", type=int, help="Lowest taxon id of the species ")
    parser.add_argument("--taxonomy_update", action="store_true", help="Update taxonomy table")

    args = parser.parse_args()

    logging.info("Loading file: %s", args.json_path)

    if not args.taxonomy_update:
        update_species_metadata(args)

    # Update taxonomy table, this is used as a standalone script
    if args.taxonomy_update and args.taxon_id:
        update_taxonomy_tables(args)


if __name__ == "__main__":
    main()
