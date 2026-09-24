# metadata_app/backend/app/api/routes/taxonomy.py

import logging
from typing import Dict, List, Optional
from fastapi import APIRouter, HTTPException

from metadata_app.backend.app.services.taxonomy_service import (
	load_clade_data,
	assign_clade_and_species,
	get_descendant_taxa
)

# Create router for taxonomy endpoints
taxonomy = APIRouter(
	prefix="/taxonomy",
	tags=["taxonomy"],
	responses={404: {"description": "Not found"}}
)


@taxonomy.get("/clade-data")
async def get_clade_data():
	"""
	Get the clade settings data.

	Returns:
		Dictionary with clade settings
	"""
	try:
		clade_data = load_clade_data()
		return {"clade_settings": clade_data}

	except Exception as e:
		logging.error(f"Error loading clade data: {str(e)}")
		raise HTTPException(status_code=500, detail=f"Internal server error: {str(e)}")


@taxonomy.get("/descendants/{taxon_id}")
async def get_descendants(taxon_id: int):
	"""
	Get all descendant taxon IDs for the given parent taxon ID.

	Args:
		taxon_id: The parent taxon ID (e.g., 40674 for Mammalia)

	Returns:
		List of descendant taxon IDs
	"""
	try:
		descendant_ids = get_descendant_taxa(taxon_id)
		return {"taxon_id": taxon_id, "descendants": list(descendant_ids)}

	except Exception as e:
		logging.error(f"Error retrieving descendant taxa: {str(e)}")
		raise HTTPException(status_code=500, detail=f"Internal server error: {str(e)}")


@taxonomy.post("/assign-clade")
async def assign_clade(request: Dict):
	"""
	Assign internal clade, species taxon ID, and genus taxon ID based on taxonomy.

	Args:
		request: Dictionary containing lowest_taxon_id and taxonomy_dict

	Returns:
		Dictionary with clade assignment information
	"""
	try:
		if "lowest_taxon_id" not in request or "taxonomy_dict" not in request:
			raise HTTPException(
				status_code=400,
				detail="Request must include 'lowest_taxon_id' and 'taxonomy_dict'"
			)

		lowest_taxon_id = request["lowest_taxon_id"]
		taxonomy_dict = request["taxonomy_dict"]

		# Optional parameters
		chordata_taxon_id = request.get("chordata_taxon_id", 7711)
		human_taxon_id = request.get("human_taxon_id", 9606)

		clade_data = load_clade_data()

		internal_clade, species_taxon_id, genus_taxon_id, pipeline = assign_clade_and_species(
			lowest_taxon_id,
			clade_data,
			taxonomy_dict,
			chordata_taxon_id,
			human_taxon_id
		)

		return {
			"lowest_taxon_id": lowest_taxon_id,
			"internal_clade": internal_clade,
			"species_taxon_id": species_taxon_id,
			"genus_taxon_id": genus_taxon_id,
			"pipeline": pipeline
		}

	except Exception as e:
		logging.error(f"Error assigning clade: {str(e)}")
		raise HTTPException(status_code=500, detail=f"Internal server error: {str(e)}")


@taxonomy.post("/batch-assign-clade")
async def batch_assign_clade(request: Dict):
	"""
	Assign clades for multiple taxon IDs in batch.

	Args:
		request: Dictionary containing taxon_ids list and taxonomy_dict

	Returns:
		Dictionary with clade assignments for each taxon ID
	"""
	try:
		if "taxon_ids" not in request or "taxonomy_dict" not in request:
			raise HTTPException(
				status_code=400,
				detail="Request must include 'taxon_ids' and 'taxonomy_dict'"
			)

		taxon_ids = request["taxon_ids"]
		taxonomy_dict = request["taxonomy_dict"]

		# Optional parameters
		chordata_taxon_id = request.get("chordata_taxon_id", 7711)
		human_taxon_id = request.get("human_taxon_id", 9606)

		clade_data = load_clade_data()

		results = {}
		for taxon_id in taxon_ids:
			internal_clade, species_taxon_id, genus_taxon_id, pipeline = assign_clade_and_species(
				taxon_id,
				clade_data,
				taxonomy_dict,
				chordata_taxon_id,
				human_taxon_id
			)

			results[taxon_id] = {
				"internal_clade": internal_clade,
				"species_taxon_id": species_taxon_id,
				"genus_taxon_id": genus_taxon_id,
				"pipeline": pipeline
			}

		return {"assignments": results}

	except Exception as e:
		logging.error(f"Error in batch clade assignment: {str(e)}")
		raise HTTPException(status_code=500, detail=f"Internal server error: {str(e)}")