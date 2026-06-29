import json
import logging
from pathlib import Path
import time
import requests


def load_clade_data():
    """Load clade settings JSON reliably."""
    base_dir = Path(__file__).resolve().parents[4]  # adjust if needed
    json_file = base_dir / "metadata_app/backend/data/clade_settings.json"

    if not json_file.exists():
        logging.error(f"Clade config not found: {json_file}")
        raise FileNotFoundError(f"{json_file} not found")

    with open(json_file, "r") as f:
        logging.info("Loading clade settings json file.")
        return json.load(f)


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

    # Build a quick mapping taxon_class -> taxon_class_id
    taxon_class_map = {
        t["taxon_class"]: t["taxon_class_id"] for t in taxonomy_hierarchy
    }

    species_taxon_id = taxon_class_map.get("species")
    genus_taxon_id = taxon_class_map.get("genus")
    species_taxon_id_int = int(species_taxon_id) if species_taxon_id is not None else None
    is_human = lowest_taxon_id == human_taxon_id or species_taxon_id_int == human_taxon_id

    # Precompute taxon_id → clade_name mapping
    clade_lookup = {
        int(details["taxon_id"]): clade_name
        for clade_name, details in clade_data.items()
        if details.get("taxon_id")
    }

    # Assign internal clade
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
        "retmax": 100000,  # Fetch in chunks
        "retstart": 0,
        "tool": "your_tool_name",
        "email": "your_email@example.com",
    }

    taxon_ids = set()

    while True:
        response = requests.get(base_url, params=params)
        if response.status_code != 200:
            logging.error(
                f"Error retrieving taxonomic data from NCBI. HTTP {response.status_code}."
            )
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
