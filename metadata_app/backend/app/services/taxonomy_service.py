import json
import logging
import time
from pathlib import Path

import requests


BACKEND_ROOT = Path(__file__).resolve().parents[2]
CLADE_DATA_FILE = BACKEND_ROOT / "data" / "clade_settings.json"


def load_clade_data():
    """Load clade settings JSON reliably."""
    if not CLADE_DATA_FILE.exists():
        logging.error("Clade config not found: %s", CLADE_DATA_FILE)
        raise FileNotFoundError(f"{CLADE_DATA_FILE} not found")

    with CLADE_DATA_FILE.open("r", encoding="utf-8") as handle:
        logging.info("Loading clade settings json file.")
        return json.load(handle)


def assign_clade_and_species(
    lowest_taxon_id,
    clade_data,
    taxonomy_dict,
    human_taxon_id=9606,
):
    """
    Assign clade, species_id, genus_id, and pipeline based on taxonomy efficiently.
    vert_taxon_id_set should be preloaded once for all records.
    """

    lowest_taxon_id = int(lowest_taxon_id)
    human_taxon_id = int(human_taxon_id)

    taxonomy_hierarchy = taxonomy_dict.get(str(lowest_taxon_id)) or taxonomy_dict.get(
        lowest_taxon_id, []
    )

    if not taxonomy_hierarchy:
        logging.warning(f"No taxonomy hierarchy found for taxon_id {lowest_taxon_id}")
        if lowest_taxon_id == human_taxon_id:
            return "human", human_taxon_id, None, "hprc"
        return "Unassigned", None, None, "anno"

    taxon_class_map = {
        taxon["taxon_class"]: taxon["taxon_class_id"]
        for taxon in taxonomy_hierarchy
    }

    species_taxon_id = taxon_class_map.get("species")
    genus_taxon_id = taxon_class_map.get("genus")
    species_taxon_id_int = (
        int(species_taxon_id) if species_taxon_id is not None else None
    )
    is_human = lowest_taxon_id == human_taxon_id or species_taxon_id_int == human_taxon_id

    clade_lookup = {
        int(details["taxon_id"]): clade_name
        for clade_name, details in clade_data.items()
        if details.get("taxon_id")
    }

    internal_clade = "Unassigned"
    for taxon_class in [
        "species",
        "genus",
        "family",
        "order",
        "class",
        "phylum",
        "kingdom",
    ]:
        taxon_id = taxon_class_map.get(taxon_class)
        if taxon_id is not None and int(taxon_id) in clade_lookup:
            internal_clade = clade_lookup[int(taxon_id)]
            break

    if is_human:
        internal_clade = "human"
        pipeline = "hprc"
    else:
        clade_settings = clade_data.get(internal_clade, {})
        pipeline = "main" if clade_settings.get("sanity_set") else "anno"

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
        "retmax": 100000,
        "retstart": 0,
        "tool": "your_tool_name",
        "email": "your_email@example.com",
    }

    taxon_ids = set()

    while True:
        response = requests.get(base_url, params=params)
        if response.status_code != 200:
            logging.error(
                "Error retrieving taxonomic data from NCBI. HTTP %s.",
                response.status_code,
            )
            break

        try:
            result = response.json()
            batch_ids = result.get("esearchresult", {}).get("idlist", [])
            logging.info("Descendant taxon ID lookup successful.")
            if not batch_ids:
                break

            taxon_ids.update(batch_ids)
            params["retstart"] += params["retmax"]
            time.sleep(0.5)

        except Exception as exc:
            logging.error("Error processing NCBI response: %s", exc)
            break

    return taxon_ids
