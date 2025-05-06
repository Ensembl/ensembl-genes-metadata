from fastapi import APIRouter, HTTPException
from metadata_app.backend.app.models.assembly_schemas import AssemblyFilterRequest
from metadata_app.backend.app.services.assembly_service import get_filtered_assemblies

router = APIRouter()

@router.post("/assemblies/filter")
def filter_assemblies(filters: AssemblyFilterRequest):
    result = get_filtered_assemblies(
        bioproject_id=filters.bioproject_id,
        metric_thresholds=filters.metric_thresholds or {},
        asm_level=filters.asm_level,
        asm_type=filters.asm_type,
        release_date=filters.release_date,
        taxon_id=filters.taxon_id,
        current=filters.current,
        pipeline=filters.pipeline,
        transc=filters.transc,
        transc_ena=filters.transc_ena,
        non_annotated = filters.non_annotated
    )

    if isinstance(result[0], str):  # Error string
        raise HTTPException(status_code=400, detail=result[0])

    df_wide, df_main, df_gca_list, taxonomy_dict = result

    return {
        "df_main": df_main.to_dict(orient="records"),
        "downloadables": {
            "df_main": df_main.to_csv(index=False),
            "df_wide": df_wide.to_csv(index=False),
            "gca_list": df_gca_list.to_csv(index=False)
        }
    }