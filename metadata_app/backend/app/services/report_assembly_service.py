import logging

import numpy as np
import pandas as pd
from fastapi import HTTPException

from metadata_app.backend.app.core.database import get_db_connection
from metadata_app.backend.app.services.assembly_service import (
    get_filtered_assemblies as get_filtered_assemblies_from_service,
)
from metadata_app.backend.app.services.get_transcriptomic_data_ENA_service import (
    add_data_from_ena,
)
from metadata_app.backend.app.services.taxonomy_service import (
    assign_clade_and_species,
    load_clade_data,
)

REPORT_ASSEMBLY_MAIN_COLUMNS = [
    "associated_project",
    "gca",
    "scientific_name",
    "release_date",
    "lowest_taxon_id",
    "genus_taxon_id",
    "transcriptomic_evidence",
    "internal_clade",
    "asm_type",
    "asm_name",
    "refseq_accession",
    "asm_level",
    "contig_n50",
    "total_sequence_length",
]


def check_dataframe_not_empty(df, description, raise_404=True):
    if df.empty:
        error_msg = f"DataFrame is empty: {description}"
        logging.error(error_msg)
        if raise_404:
            raise HTTPException(status_code=404, detail=f"No data found: {description}")
        return False

    logging.info(f"DataFrame validated - {description}: {len(df)} rows")
    return True


def _join_unique(values):
    return ", ".join(
        sorted({str(value) for value in values if pd.notna(value) and value != ""})
    )


def _add_transcriptomic_evidence(df_wide):
    if "transcriptomic_evidence" in df_wide.columns:
        return df_wide

    df_wide = df_wide.copy()
    lowest_col = "short_read_paired_end_illumina_lowest"
    genus_col = "short_read_paired_end_illumina"

    if lowest_col in df_wide.columns and genus_col in df_wide.columns:
        lowest_reads = pd.to_numeric(df_wide[lowest_col], errors="coerce")
        genus_reads = pd.to_numeric(df_wide[genus_col], errors="coerce")
        df_wide["transcriptomic_evidence"] = np.where(
            (lowest_reads != 0) | ((lowest_reads == 0) & (genus_reads >= 2)),
            "yes",
            "no",
        )
    else:
        df_wide["transcriptomic_evidence"] = "not checked"

    return df_wide


def _apply_report_filters(df_wide, candidate, start_date):
    df_wide = df_wide.copy()

    if start_date:
        df_wide["release_date"] = pd.to_datetime(
            df_wide["release_date"], errors="coerce"
        )
        df_wide = df_wide[df_wide["release_date"] <= pd.to_datetime(start_date)]
        check_dataframe_not_empty(df_wide, "assemblies before report start date")

    if candidate:
        if "contig_n50" not in df_wide.columns:
            raise HTTPException(
                status_code=404,
                detail="No contig_n50 metric found for candidate assembly filtering.",
            )
        df_wide["contig_n50"] = pd.to_numeric(df_wide["contig_n50"], errors="coerce")
        df_wide = df_wide[
            (~df_wide["asm_level"].isin(["Contig", "Scaffold"]))
            & (df_wide["contig_n50"] >= 100000)
        ]
        logging.info(
            "Filtered for annotation candidates (n50>=100000 and asm_level!='Contig')"
        )
        check_dataframe_not_empty(
            df_wide, "candidate assemblies (N50 >= 100000, non-contig)"
        )

    return df_wide


def _build_report_tables(df_wide):
    df_wide = _add_transcriptomic_evidence(df_wide)

    missing_columns = [
        column for column in REPORT_ASSEMBLY_MAIN_COLUMNS if column not in df_wide
    ]
    if missing_columns:
        raise HTTPException(
            status_code=500,
            detail=f"Assembly report is missing expected columns: {', '.join(missing_columns)}",
        )

    agg_dict = {
        "bioproject_id": _join_unique,
        "associated_project": _join_unique,
    }

    for col in df_wide.columns:
        if col not in ["gca", "bioproject_id", "associated_project"]:
            agg_dict[col] = "first"

    rep_asm_wide = df_wide.groupby("gca", as_index=False).agg(agg_dict)
    check_dataframe_not_empty(
        rep_asm_wide, "wide results table after deduplication and cleanup"
    )

    rep_asm_main = df_wide[REPORT_ASSEMBLY_MAIN_COLUMNS]
    check_dataframe_not_empty(rep_asm_main, "main results table after deduplication")
    logging.info("Created main table")

    return rep_asm_wide, rep_asm_main


def generate_tables(
    bioproject_id,
    candidate,
    taxon_id,
    transc,
    transc_ena,
    start_date,
    end_date,
    group_name,
    non_annotated=True,
):
    logging.info(
        f"Generating tables for end date: {end_date}, start date: {start_date}, group name: {group_name}, taxon id: {taxon_id}, bioproject id: {bioproject_id}, non_annotated: {non_annotated}"
    )

    try:
        assembly_result = get_filtered_assemblies_from_service(
            bioproject_id=bioproject_id,
            metric_thresholds={},
            asm_level=None,
            asm_type=None,
            release_date=end_date,
            taxon_id=taxon_id,
            current=True,
            transc=transc,
            transc_ena=transc_ena,
            non_annotated=non_annotated,
            group_name=group_name,
            gca=None,
        )
        if isinstance(assembly_result[0], str):
            return assembly_result

        df_wide, _, _ = assembly_result
        df_wide = _apply_report_filters(df_wide, candidate, start_date)
        rep_asm_wide, rep_asm_main = _build_report_tables(df_wide)
    except HTTPException:
        logging.error("HTTPException raised during assembly filtering")
        raise
    except Exception:
        logging.error(
            "Unexpected error occurred during assembly filtering", exc_info=True
        )
        raise HTTPException(status_code=500, detail="An unexpected error occurred.")

    check_dataframe_not_empty(rep_asm_wide, "wide assembly results from filtering")
    check_dataframe_not_empty(rep_asm_main, "main assembly results from filtering")

    df_gca_list = rep_asm_wide[["gca"]]
    check_dataframe_not_empty(df_gca_list, "GCA list")
    logging.info("Created gca_list")

    project_report = (
        rep_asm_wide[["gca", "associated_project"]]
        .groupby("associated_project")
        .size()
        .reset_index(name="count")
    )
    check_dataframe_not_empty(project_report, "project report summary")

    transc_data = (
        rep_asm_wide[["gca", "transcriptomic_evidence"]]
        .groupby("transcriptomic_evidence")
        .size()
        .reset_index(name="count")
    )
    check_dataframe_not_empty(transc_data, "transcriptomic data summary")

    num_unique_taxa = rep_asm_wide["lowest_taxon_id"].nunique()
    if num_unique_taxa == 0:
        logging.error("No unique taxa found")
        raise HTTPException(
            status_code=404, detail="No unique taxa found in the results"
        )

    top_3_taxa = (
        rep_asm_wide.groupby(["scientific_name"])
        .size()
        .reset_index(name="count")
        .sort_values(by="count", ascending=False)
        .head(3)
    )
    check_dataframe_not_empty(top_3_taxa, "top 3 taxa summary")

    asm_type_group = (
        rep_asm_wide[["gca", "asm_type"]]
        .groupby("asm_type")
        .size()
        .reset_index(name="count")
    )

    asm_level_group = (
        rep_asm_wide[["gca", "asm_level"]]
        .groupby("asm_level")
        .size()
        .reset_index(name="count")
    )
    check_dataframe_not_empty(asm_type_group, "assembly type summary")

    clade_group = (
        rep_asm_wide[["gca", "internal_clade"]]
        .groupby("internal_clade")
        .size()
        .reset_index(name="count")
    )
    check_dataframe_not_empty(clade_group, "clade group summary")

    asm_length = rep_asm_wide[["gca", "total_sequence_length"]].copy()
    check_dataframe_not_empty(asm_length, "assembly length data")

    asm_length = asm_length.sort_values(
        by="total_sequence_length", ascending=False
    ).reset_index(drop=True)
    asm_length["total_sequence_length"] = pd.to_numeric(
        asm_length["total_sequence_length"], errors="coerce"
    )
    asm_length["total_sequence_length_Gb"] = asm_length["total_sequence_length"] / 1e9
    asm_length = asm_length[["gca", "total_sequence_length_Gb"]]
    check_dataframe_not_empty(asm_length, "processed assembly length data")
    logging.info(f"Assembly length data processed successfully: {len(asm_length)} rows")

    transc_cols = [
        col for col in rep_asm_wide.columns if col.endswith("_transc_assess_date")
    ]
    if not transc_cols:
        transc_reg_count = "not checked"
    else:
        total_rows = len(rep_asm_wide)
        rows_with_transc = (
            rep_asm_wide[transc_cols]
            .apply(
                lambda row: any(
                    pd.notna(val) and str(val).strip() != "" for val in row
                ),
                axis=1,
            )
            .sum()
        )
        transc_reg_count = f"{(rows_with_transc / total_rows) * 100:.1f}%"

    logging.info("Transforming out of range float values that are not JSON compliant")
    rep_asm_wide = rep_asm_wide.apply(
        lambda col: col.fillna("") if col.dtype == "object" else col
    )
    project_report = project_report.apply(
        lambda col: col.fillna("") if col.dtype == "object" else col
    )
    rep_asm_main = rep_asm_main.apply(
        lambda col: col.fillna("") if col.dtype == "object" else col
    )
    top_3_taxa = top_3_taxa.apply(
        lambda col: col.fillna("") if col.dtype == "object" else col
    )
    asm_type_group = asm_type_group.apply(
        lambda col: col.fillna("") if col.dtype == "object" else col
    )
    asm_level_group = asm_level_group.apply(
        lambda col: col.fillna("") if col.dtype == "object" else col
    )
    clade_group = clade_group.apply(
        lambda col: col.fillna("") if col.dtype == "object" else col
    )
    asm_length = asm_length.apply(
        lambda col: col.fillna("") if col.dtype == "object" else col
    )
    transc_data = transc_data.apply(
        lambda col: col.fillna("") if col.dtype == "object" else col
    )

    check_dataframe_not_empty(project_report, "final project report")
    check_dataframe_not_empty(top_3_taxa, "final top 3 taxa")
    check_dataframe_not_empty(asm_type_group, "final assembly type group")
    check_dataframe_not_empty(asm_level_group, "final assembly level group")
    check_dataframe_not_empty(clade_group, "final clade group")
    check_dataframe_not_empty(asm_length, "final assembly length")
    check_dataframe_not_empty(transc_data, "final transcriptomic data")
    check_dataframe_not_empty(df_gca_list, "final GCA list")
    check_dataframe_not_empty(rep_asm_wide, "final wide assembly results")
    check_dataframe_not_empty(rep_asm_main, "final main assembly results")

    logging.info("All tables generated successfully")

    return (
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
    )


def generate_overview():
    try:
        with get_db_connection("meta") as conn:
            cursor = conn.cursor()

            query = """
                SELECT 
                    COALESCE(g.group_name, mb.bioproject_name) AS project_name,
                    a.assembly_id,
                    a.asm_level,
                    m.metrics_name,
                    m.metrics_value,
                    a.lowest_taxon_id,
                    s.species_taxon_id,
                    a.is_current,
                    t.taxon_class_id AS genus_taxon_id,
                    gb.gb_status AS genebuild_status
                FROM assembly a
                LEFT JOIN bioproject b ON a.assembly_id = b.assembly_id
                LEFT JOIN main_bioproject mb ON mb.bioproject_id = b.bioproject_id
                LEFT JOIN assembly_metrics m ON b.assembly_id = m.assembly_id
                JOIN taxonomy t ON a.lowest_taxon_id = t.lowest_taxon_id
                JOIN species s ON a.lowest_taxon_id = s.lowest_taxon_id
                LEFT JOIN custom_group g
                    ON (
                        (g.group_type = 'taxon' AND a.lowest_taxon_id = g.item)
                        OR
                        (g.group_type = 'assembly' AND a.gca_chain = g.item)
                    )
                LEFT JOIN genebuild_status gb ON a.assembly_id = gb.assembly_id
                WHERE (mb.bioproject_id IS NOT NULL OR g.group_name IS NOT NULL)
                  AND t.taxon_class = "genus"
                  AND m.metrics_name IN ('contig_n50');
            """

            cursor.execute(query)
            results = cursor.fetchall()

        if not results:
            raise HTTPException(
                status_code=404, detail="No assemblies found matching criteria."
            )

        df = pd.DataFrame(results)
        logging.info(f"Fetched {len(df)} rows from database")
        print(f"Projects before pivot: {df['project_name'].nunique()}")

        df["genebuild_status"] = df["genebuild_status"].fillna("not_annotated")
        df = df.drop_duplicates(
            subset=["project_name", "assembly_id", "lowest_taxon_id"], keep="first"
        )

        lowest_taxon_ids = {
            row["lowest_taxon_id"]
            for row in results
            if row.get("lowest_taxon_id") is not None
        }
        taxonomy_dict = {}
        if lowest_taxon_ids:
            with get_db_connection("meta") as conn:
                cursor = conn.cursor()
                taxonomy_query = """
                    SELECT lowest_taxon_id, taxon_class_id, taxon_class
                    FROM taxonomy
                    WHERE lowest_taxon_id IN ({})
                """.format(",".join(["%s"] * len(lowest_taxon_ids)))
                cursor.execute(taxonomy_query, tuple(lowest_taxon_ids))
                taxonomy_results = cursor.fetchall()

            for row in taxonomy_results:
                lowest_taxon_id = row["lowest_taxon_id"]
                taxonomy_dict.setdefault(lowest_taxon_id, []).append(
                    {
                        "taxon_class_id": row["taxon_class_id"],
                        "taxon_class": row["taxon_class"],
                    }
                )

        clade_data = load_clade_data()
        df[["internal_clade", "species_taxon_id", "genus_taxon_id", "pipeline"]] = (
            df["lowest_taxon_id"].apply(
                lambda x: pd.Series(
                    assign_clade_and_species(x, clade_data, taxonomy_dict)
                )
            )
        )

        df_wide = df.pivot_table(
            index=[
                "project_name",
                "assembly_id",
                "asm_level",
                "is_current",
                "lowest_taxon_id",
                "species_taxon_id",
                "genus_taxon_id",
                "genebuild_status",
                "pipeline",
            ],
            columns="metrics_name",
            values="metrics_value",
            aggfunc="first",
            fill_value=np.nan,
        ).reset_index()

        print(f"Projects after pivot: {df_wide['project_name'].nunique()}")

        transcriptomic_df = add_data_from_ena(df_wide)
        if transcriptomic_df is not None and not transcriptomic_df.empty:
            df_wide = df_wide.merge(
                transcriptomic_df[["taxon_id", "short_read_paired_end_illumina"]],
                left_on="genus_taxon_id",
                right_on="taxon_id",
                how="left",
            )
            df_wide["transcriptomic_evidence"] = np.where(
                df_wide["short_read_paired_end_illumina"].fillna(0) > 2, "yes", "no"
            )
        else:
            df_wide["transcriptomic_evidence"] = "no"

        df_wide["is_annotation_candidate"] = (
            df_wide["asm_level"].isin(["Complete genome", "Chromosome"])
            & (df_wide["contig_n50"].astype(float) > 100000)
            & (df_wide["transcriptomic_evidence"] == "yes")
            & (df_wide["is_current"] == "current")
        )

        df_wide["unannotated"] = (
            df_wide["asm_level"].isin(["Complete genome", "Chromosome"])
            & (df_wide["contig_n50"].astype(float) > 100000)
            & (df_wide["transcriptomic_evidence"] == "yes")
            & (df_wide["is_current"] == "current")
            & (df_wide["genebuild_status"] == "not_annotated")
        )

        df_wide["unannotated_main"] = (
            df_wide["unannotated"] & (df_wide["pipeline"] == "main")
        )
        df_wide["unannotated_anno"] = (
            df_wide["unannotated"] & (df_wide["pipeline"] == "anno")
        )

        summary = df_wide.groupby("project_name", as_index=False).agg(
            total_assemblies=("assembly_id", "nunique"),
            annotation_candidates=("is_annotation_candidate", "sum"),
            unannotated=("unannotated", "sum"),
            unannotated_main=("unannotated_main", "sum"),
            unannotated_anno=("unannotated_anno", "sum"),
            in_progress=(
                "genebuild_status",
                lambda x: x.isin(["in_progress", "complete", "pre_released"]).sum(),
            ),
            live=("genebuild_status", lambda x: (x == "live").sum()),
        )

        logging.info(f"Generated overview summary for {len(summary)} projects")
        return summary.to_dict(orient="records")

    except Exception as e:
        logging.exception("Error generating overview table")
        raise HTTPException(status_code=500, detail=str(e))


def project_per_year():
    try:
        with get_db_connection("meta") as conn:
            cursor = conn.cursor()

            query = """
                SELECT 
                    COALESCE(g.group_name, mb.bioproject_name) AS project_name,
                    a.assembly_id,
                    a.release_date,
                    gb.release_date AS release_date_anno,
                    gb.gb_status, 
                    a.is_current
                FROM assembly a
                LEFT JOIN bioproject b ON b.assembly_id = a.assembly_id
                LEFT JOIN main_bioproject mb ON mb.bioproject_id = b.bioproject_id
                LEFT JOIN custom_group g
                    ON (
                        (g.group_type = 'taxon' AND a.lowest_taxon_id = g.item)
                        OR
                        (g.group_type = 'assembly' AND a.gca_chain = g.item)
                    )
                LEFT JOIN genebuild_status gb ON a.assembly_id = gb.assembly_id
                WHERE (mb.bioproject_id IS NOT NULL OR g.group_name IS NOT NULL)
            """

            cursor.execute(query)
            results = cursor.fetchall()

        if not results:
            raise HTTPException(
                status_code=404, detail="No assemblies found matching criteria."
            )

        df = pd.DataFrame(results)
        logging.info(f"Fetched {len(df)} rows from database")
        print(f"Projects per year before: {df['project_name'].nunique()}")

        df["gb_status"] = df["gb_status"].fillna("not_annotated")
        df = df.drop_duplicates(subset=["project_name", "assembly_id"], keep="first")
        print(
            f"Projects per year after drop duplicates: {df['project_name'].nunique()}"
        )

        df["release_date"] = pd.to_datetime(df["release_date"], errors="coerce")
        df["release_date_anno"] = pd.to_datetime(
            df["release_date_anno"], errors="coerce"
        )
        print(f"Projects per year after datetime: {df['project_name'].nunique()}")

        df_assembly = df
        df_assembly["release_year"] = df_assembly["release_date"].dt.year
        df_assembly = df_assembly[df_assembly["release_year"] >= 2019]

        df_live = df
        df_live["release_year_anno"] = df_live["release_date_anno"].dt.year

        print(
            f"Projects per year after extract year: {df_live['project_name'].nunique()}"
        )
        print(df_live)

        df_assembly = df_assembly[df_assembly["is_current"] == "current"]
        df_assembly = df_assembly.groupby(
            ["project_name", "release_year"], as_index=False
        ).agg(total_assemblies=("assembly_id", "nunique"))

        df_live = df_live[df_live["gb_status"] == "live"].copy()
        print(
            f"Projects per year after live filter by annotations: {df_live['project_name'].nunique()}"
        )

        df_annotations = (
            df_live.groupby(["project_name", "release_year_anno"], as_index=False)
            .agg(live_annotations=("assembly_id", "nunique"))
            .rename(columns={"release_year_anno": "release_year"})
        )
        print(
            f"Projects per year after group by annotations: {df_annotations['project_name'].nunique()}"
        )

        df_assembly = df_assembly.pivot_table(
            index="release_year",
            columns="project_name",
            values="total_assemblies",
            fill_value=0,
        ).reset_index()
        df_assembly.columns.name = None

        df_annotations = df_annotations.pivot_table(
            index="release_year",
            columns="project_name",
            values="live_annotations",
            fill_value=0,
        ).reset_index()
        df_annotations.columns.name = None

        df_assembly["release_year"] = df_assembly["release_year"].astype(int)
        df_annotations["release_year"] = df_annotations["release_year"].astype(int)

        df_assembly = df_assembly.to_dict(orient="records")
        df_annotations = df_annotations.to_dict(orient="records")
        print("Assemblies bar")
        print(df_assembly)
        print("Annotations bar")
        print(df_annotations)

        return df_assembly, df_annotations

    except Exception as e:
        logging.exception("Error generating project per year tables")
        raise HTTPException(status_code=500, detail=str(e))
