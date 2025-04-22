from fastapi import FastAPI, Query
from pydantic import BaseModel
from typing import List, Optional
import logging
import
app = FastAPI()

# Pydantic model to validate incoming request data
class FilterRequest(BaseModel):
    bioproject_id: List[str]
    gc_percent: Optional[float] = None
    total_sequence_length: Optional[float] = None
    contig_n50: Optional[float] = None
    number_of_contigs: Optional[float] = None
    number_of_scaffolds: Optional[float] = None
    scaffold_n50: Optional[float] = None
    genome_coverage: Optional[float] = None
    release_date: Optional[str] = None
    asm_level: Optional[List[str]] = None
    asm_type: Optional[List[str]] = None
    taxon_id: Optional[int] = None
    current: Optional[int] = 0
    pipeline: Optional[List[str]] = None
    output_dir: str = "./"

@app.get("/")
def healthcheck():
    return {"status": "ok"}

@app.post("/filter_assemblies/")
async def filter_assemblies(request: FilterRequest):
    try:
        # Call the function to process the request
        result = filter_assemblies_and_check_transcriptomic_data(
            bioproject_id=request.bioproject_id,
            metric_thresholds={key: value for key, value in request.dict().items() if value is not None},
            asm_level=request.asm_level,
            asm_type=request.asm_type,
            release_date=request.release_date,
            taxon_id=request.taxon_id,
            current=request.current,
            pipeline=request.pipeline,
            output_dir=request.output_dir
        )
        return result
    except Exception as e:
        logging.error(f"Error processing request: {str(e)}")
        raise HTTPException(status_code=500, detail="Internal server error")