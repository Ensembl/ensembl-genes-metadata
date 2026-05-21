import logging

from fastapi import APIRouter, HTTPException
from metadata_app.backend.app.models.report_asm_schemas import ReportFilterRequest
from metadata_app.backend.app.services.report_assembly_service import (
    generate_tables,
    generate_overview,
    project_per_year,
)

report = APIRouter()


@report.post("/report/asm/filter")
def filter_assemblies(filters: ReportFilterRequest):
    result = generate_tables(
        bioproject_id=filters.bioproject_id,
        start_date=filters.start_date,
        end_date=filters.end_date,
        group_name=filters.group_name,
        taxon_id=filters.taxon_id,
        candidate=filters.candidate,
        transc=filters.transc,
        transc_ena=filters.transc_ena,
        non_annotated=(
            filters.non_annotated if filters.non_annotated is not None else True
        ),
    )

    if isinstance(result[0], str):  # Error string
        raise HTTPException(status_code=400, detail=result[0])

    (
        project_report,
        num_unique_taxa,
        transc_reg_count,
        top_3_taxa,
        asm_type_group,
        asm_level_group,
        clade_group,
        asm_length,
        transc_data,
        df_gca_list,
        rep_asm_wide,
        rep_asm_main,
    ) = result

    return {
        "project_report": project_report.to_dict(orient="records"),
        "num_unique_taxa": {"value": num_unique_taxa},
        "transc_data": transc_data.to_dict(orient="records"),
        "top_3_taxa": top_3_taxa.to_dict(orient="records"),
        "asm_type_group": asm_type_group.to_dict(orient="records"),
        "asm_level_group": asm_level_group.to_dict(orient="records"),
        "clade_group": clade_group.to_dict(orient="records"),
        "asm_length": asm_length.to_dict(orient="records"),
        "transc_reg_count": {"value": transc_reg_count},
        "rep_asm_main": rep_asm_main.to_dict(orient="records"),
        "downloadables_report": {
            "rep_asm_main": rep_asm_main.to_csv(index=False),
            "rep_asm_wide": rep_asm_wide.to_csv(index=False),
            "gca_list": df_gca_list.to_csv(index=False),
        },
    }


@report.get("/report/asm/main")
def generate_main():
    return generate_overview()


@report.get("/report/asm/bar")
def generate_bar():
    df_assembly, df_annotations = project_per_year()
    return {"df_assembly": df_assembly, "df_annotations": df_annotations}
