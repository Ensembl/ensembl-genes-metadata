import json
import logging
import time
import requests


def load_clade_data():
    """Hardcoded path for clade settings."""
    json_file = "data/clade_settings.json"
    with open(json_file, "r") as f:
        logging.info("Loading clade settings json file.")
        return json.load(f)

def read_ids_as_set(file_path: str) -> set:
    """
    Reads a text file with one ID per line and returns a set of integers.

    :param file_path: Path to the text file
    :return: A set of unique IDs as integers
    """
    with open(file_path, "r") as f:
        return {int(line.strip()) for line in f if line.strip()}


def assign_clade_and_species(lowest_taxon_id, clade_data, taxonomy_dict, human_taxon_id=9606):
    """Assign internal clade and species taxon ID based on taxonomy using the provided clade data,
       and check if the taxon ID is a descendant of the vertebrata taxon ID (7742)."""

    # Convert IDs to integers to ensure consistent comparison
    lowest_taxon_id = int(lowest_taxon_id)
    human_taxon_id = int(human_taxon_id)

    # Add debug logging to see what's being passed
    logging.debug(f"Processing taxon ID: {lowest_taxon_id}")
    logging.debug(f"Taxonomy dict keys available: {list(taxonomy_dict.keys())}")

    # Retrieve the taxonomy hierarchy from the passed dictionary
    taxonomy_hierarchy = taxonomy_dict.get(str(lowest_taxon_id)) or taxonomy_dict.get(lowest_taxon_id, [])

    logging.debug(f"Taxonomy hierarchy for {lowest_taxon_id}: {taxonomy_hierarchy}")

    if not taxonomy_hierarchy:
        logging.warning(f"Taxonomy hierarchy not found for taxon ID {lowest_taxon_id}")
        return "Unassigned", None, None, None

    species_taxon_id = None
    genus_taxon_id = None
    pipeline = None
    internal_clade = "Unassigned"  # Default value if no clade is found

    # Define the taxon classes in hierarchical order
    taxon_classes_order = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']

    # First pass: Set species and genus taxon IDs if available
    for taxon in taxonomy_hierarchy:
        if taxon['taxon_class'] == 'species':
            species_taxon_id = taxon['taxon_class_id']
        elif taxon['taxon_class'] == 'genus':
            genus_taxon_id = taxon['taxon_class_id']

    # Second pass: Try to assign clade
    for taxon_class in taxon_classes_order:
        matching_taxon = next((t for t in taxonomy_hierarchy if t['taxon_class'] == taxon_class), None)

        if matching_taxon:
            taxon_class_id = matching_taxon['taxon_class_id']

            # Log clade data for debugging
            logging.debug(f"Checking taxon class: {taxon_class}, id: {taxon_class_id}")
            logging.debug(f"Available clades: {list(clade_data.keys())}")

            # Check for matching taxon_id in clade settings
            for clade_name, details in clade_data.items():
                clade_taxon_id = details.get("taxon_id")
                logging.debug(f"Comparing with clade {clade_name}, taxon_id: {clade_taxon_id}")

                if clade_taxon_id and int(clade_taxon_id) == int(taxon_class_id):
                    internal_clade = clade_name
                    logging.info(f"Assigned clade {clade_name} for taxon {lowest_taxon_id}")
                    break

            if internal_clade != "Unassigned":
                break

    # Check if vertebrate and assign pipeline
    vert_taxon_id_set = read_ids_as_set("data/vertebrata_taxids.txt")
    if lowest_taxon_id == human_taxon_id:
        pipeline = "hprc"
    elif lowest_taxon_id in vert_taxon_id_set:
        pipeline = "main"
    else:
        pipeline = "anno"

    # Log the assignment results for debugging
    logging.info(
        f"Taxon {lowest_taxon_id}: clade={internal_clade}, species_id={species_taxon_id}, genus_id={genus_taxon_id}, pipeline={pipeline}")

    return internal_clade, species_taxon_id, genus_taxon_id, pipeline

def get_descendant_taxa(taxon_id):
    """
    Retrieves all descendant taxon IDs under the given taxon ID using NCBI E-utilities with pagination.

    :param taxon_id: The parent taxon ID (e.g., 40674 for Mammalia)
    :return: A set of descendant taxon IDs
    """
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "taxonomy",
        "term": f"txid{taxon_id}[Subtree]",
        "retmode": "json",
        "retmax": 100000,  # Fetch in chunks
        "retstart": 0,
        "tool": "your_tool_name",
        "email": "your_email@example.com"
    }

    taxon_ids = set()

    while True:
        response = requests.get(base_url, params=params)
        if response.status_code != 200:
            logging.error(f"Error retrieving taxonomic data from NCBI. HTTP {response.status_code}.")
            break

        try:
            result = response.json()
            batch_ids = result.get("esearchresult", {}).get("idlist", [])
            logging.info(f"Descendant taxon ID lookup successful.")
            if not batch_ids:
                break  # No more results

            taxon_ids.update(batch_ids)

            # Update retstart to fetch the next batch
            params["retstart"] += params["retmax"]

            # Respect NCBI rate limits
            time.sleep(0.5)  # Avoid overloading NCBI servers

        except Exception as e:
            print(f"Error processing response: {e}")
            logging.error(f"Error processing NCBI response: {e}")
            break

    return taxon_ids