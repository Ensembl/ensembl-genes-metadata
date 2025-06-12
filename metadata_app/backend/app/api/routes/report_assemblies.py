import logging

from fastapi import APIRouter, HTTPException
from metadata_app.backend.app.models.report_asm_schemas import ReportFilterRequest
from metadata_app.backend.app.services.report_assembly_service import generate_tables

report = APIRouter()

@report.post("/report/asm/filter")
def filter_assemblies(filters: ReportFilterRequest):
    result = generate_tables(
        bioproject_id=filters.bioproject_id,
        start_date=filters.start_date,
        end_date=filters.end_date,
        group_name=filters.group_name,
        taxon_id=filters.taxon_id,
        pipeline=filters.pipeline,
        candidate = filters.candidate,
        transc=filters.transc,
        transc_ena=filters.transc_ena,
    )

    if isinstance(result[0], str):  # Error string
        raise HTTPException(status_code=400, detail=result[0])

    project_report, num_unique_taxa, transc_reg_count, top_3_taxa,asm_type_group, clade_group, asm_length, transc_data, df_gca_list, rep_asm_wide, rep_asm_main = result

    return {
        "project_report": project_report.to_dict(orient="records"),
        "num_unique_taxa": {"value": num_unique_taxa},
        "transc_data": transc_data.to_dict(orient="records"),
        "top_3_taxa": top_3_taxa.to_dict(orient="records"),
        "asm_type_group": asm_type_group.to_dict(orient="records"),
        "clade_group": clade_group.to_dict(orient="records"),
        "asm_length": asm_length.to_dict(orient="records"),
        "transc_reg_count": {"value": transc_reg_count},
        "rep_asm_main": rep_asm_main.to_dict(orient="records"),


        "downloadables_report": {
            "rep_asm_main": rep_asm_main.to_csv(index=False),
            "rep_asm_wide": rep_asm_wide.to_csv(index=False),
            "gca_list": df_gca_list.to_csv(index=False)
        }
    }