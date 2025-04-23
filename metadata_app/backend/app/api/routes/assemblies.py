# app/api/routes/assemblies.py
from fastapi import APIRouter, Query, Path, HTTPException, Depends
from typing import List, Optional, Dict, Any
import pandas as pd
from metadata_app.backend.app.services.assembly_service import (
	get_filtered_assemblies,
	is_reference_genome
)

router = APIRouter()


@router.get("/filter", response_model=Dict[str, Any])
async def filter_assemblies(
		bioproject_id: Optional[List[str]] = Query(None, description="One or more BioProject IDs"),
		gc_percent: Optional[float] = Query(None, description="Min GC percent threshold"),
		total_sequence_length: Optional[float] = Query(None, description="Min total sequence length in bp"),
		contig_n50: Optional[float] = Query(None, description="Min contig N50"),
		number_of_contigs: Optional[float] = Query(None, description="Max number of contigs"),
		number_of_scaffolds: Optional[float] = Query(None, description="Max number of scaffolds"),
		scaffold_n50: Optional[float] = Query(None, description="Min scaffold N50"),
		genome_coverage: Optional[float] = Query(None, description="Min genome coverage"),
		release_date: Optional[str] = Query(None,
		                                    description="Filter assemblies released after this date (YYYY-MM-DD)"),
		asm_level: Optional[List[str]] = Query(None,
		                                       description="Assembly levels (Contig, Scaffold, Chromosome, Complete genome)"),
		asm_type: Optional[List[str]] = Query(None,
		                                      description="Assembly types (haploid, alternate-pseudohaplotype, etc.)"),
		taxon_id: Optional[int] = Query(None, description="NCBI Taxon ID to filter by"),
		reference: bool = Query(False, description="Check if GCA is a reference genome"),
		transcriptomic: bool = Query(False, description="Check for transcriptomic data"),
		current: bool = Query(False, description="Filter for current assemblies only"),
		pipeline: Optional[List[str]] = Query(None, description="Pipeline to filter by (anno, main, hprc)")
):
	"""
	Filter genomic assemblies based on various criteria.

	Returns a dictionary containing:
	- assembly_data: Filtered assembly metrics
	- summary: Summary statistics of the filtered data
	- info: Additional assembly information
	- gca_list: List of GCA accessions
	"""
	try:
		# Format parameters for the service function
		metric_thresholds = {}
		if gc_percent is not None:
			metric_thresholds["gc_percent"] = gc_percent
		if total_sequence_length is not None:
			metric_thresholds["total_sequence_length"] = total_sequence_length
		if contig_n50 is not None:
			metric_thresholds["contig_n50"] = contig_n50
		if number_of_contigs is not None:
			metric_thresholds["number_of_contigs"] = number_of_contigs
		if number_of_scaffolds is not None:
			metric_thresholds["number_of_scaffolds"] = number_of_scaffolds
		if scaffold_n50 is not None:
			metric_thresholds["scaffold_n50"] = scaffold_n50
		if genome_coverage is not None:
			metric_thresholds["genome_coverage"] = genome_coverage

		all_metrics = ["gc_percent", "total_sequence_length", "contig_n50",
		               "number_of_contigs", "number_of_scaffolds", "scaffold_n50",
		               "genome_coverage"]

		# Get filtered assemblies
		df, summary_df, info_result, df_gca_list, taxonomy_dict = get_filtered_assemblies(
			bioproject_id=bioproject_id,
			metric_thresholds=metric_thresholds,
			all_metrics=all_metrics,
			asm_level=asm_level,
			asm_type=asm_type,
			release_date=release_date,
			taxon_id=taxon_id,
			current=current,
			pipeline=pipeline
		)

		# Handle string result (error case)
		if isinstance(df, str):
			raise HTTPException(status_code=404, detail=df)

		# Check reference genome if requested
		if reference and not info_result.empty:
			info_result["reference_genome"] = info_result["gca"].apply(is_reference_genome)

		# Check transcriptomic data if requested
		if transcriptomic and not info_result.empty:
			# This would call the transcriptomic service
			# For now, returning without additional processing
			pass

		# Convert pandas DataFrames to dictionaries
		return {
			"assembly_data": df.to_dict(orient="records") if not isinstance(df, str) else [],
			"summary": summary_df.to_dict() if not isinstance(summary_df, str) else {},
			"info": info_result.to_dict(orient="records") if not isinstance(info_result, str) else [],
			"gca_list": df_gca_list["gca"].tolist() if not isinstance(df_gca_list, str) else []
		}
	except Exception as e:
		raise HTTPException(status_code=500, detail=f"An error occurred: {str(e)}")


@router.get("/{gca}", response_model=Dict[str, Any])
async def get_assembly_details(
		gca: str = Path(..., description="GCA accession number")
):
	"""Get detailed information for a specific GCA accession"""
	try:
		# This would call a service function to get details for a specific GCA
		# For demonstration, we'll return a placeholder
		return {
			"gca": gca,
			"message": "This endpoint would return detailed information for the specified GCA."
		}
	except Exception as e:
		raise HTTPException(status_code=500, detail=f"An error occurred: {str(e)}")