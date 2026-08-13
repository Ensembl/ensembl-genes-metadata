# app/services/home_page_service.py
import logging
from datetime import date, datetime

import numpy as np
import pandas as pd
from metadata_app.backend.app.core.database import get_db_connection
from metadata_app.backend.app.services.taxonomy_service import (
    assign_clade_and_species,
    load_clade_data,
)
import pymysql.cursors


def _json_safe(value):
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_json_safe(item) for item in value]
    if isinstance(value, tuple):
        return tuple(_json_safe(item) for item in value)
    if value is None:
        return None
    if isinstance(value, pd.Timestamp):
        return None if pd.isna(value) else value.strftime("%Y-%m-%d")
    if isinstance(value, (datetime, date)):
        return value.isoformat()
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, float):
        return value if np.isfinite(value) else None
    if isinstance(value, np.bool_):
        return bool(value)
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    return value


def _json_safe_records(frame):
    return _json_safe(frame.to_dict(orient="records"))


def _empty_handover_result():
    return ([], 0, 0, 0, [], [], [], [])


def _build_core_name(scientific_name, gca):
    if not scientific_name or not gca:
        return None

    name_parts = str(scientific_name).strip().lower().split()
    if not name_parts:
        return None

    species_part = "_".join(name_parts[:2])
    gca_part = str(gca).strip().lower().replace(".", "v").replace("_", "")
    return f"{species_part}_{gca_part}_core_114_1"


def get_ready_to_ho(genebuilder):
    try:
        with get_db_connection("meta") as conn:
            cursor = conn.cursor(pymysql.cursors.DictCursor)

            query = """
				SELECT 
					g.genebuild_status_id,
					g.assembly_id,
					g.gb_status,
					g.date_status_update,
					m.bioproject_name,
					s.scientific_name,
					a.lowest_taxon_id,
					g.annotation_method,
					CONCAT(a.gca_chain, ".", a.gca_version) AS gca,
					(
						SELECT COUNT(*)
						FROM genebuild_status g2
						WHERE g2.assembly_id = g.assembly_id
						  AND g2.genebuild_status_id <> g.genebuild_status_id
						  AND g2.gb_status <> 'abandoned'
					) AS other_annotation_count
				FROM genebuild_status g
				JOIN assembly a ON a.assembly_id = g.assembly_id
				LEFT JOIN bioproject b ON b.assembly_id = a.assembly_id
				LEFT JOIN main_bioproject m ON m.bioproject_id = b.bioproject_id
				LEFT JOIN species s ON s.lowest_taxon_id = a.lowest_taxon_id
				WHERE g.genebuilder = %s
			"""

            cursor.execute(query, (genebuilder,))
            result = cursor.fetchall()

            lowest_taxon_ids = {
                row["lowest_taxon_id"]
                for row in result
                if row.get("lowest_taxon_id") is not None
            }
            taxonomy_dict = {}
            if lowest_taxon_ids:
                taxonomy_query = """
                    SELECT lowest_taxon_id, taxon_class_id, taxon_class
                    FROM taxonomy
                    WHERE lowest_taxon_id IN ({})
                    ORDER BY FIELD(taxon_class, 'species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom');
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

        if not result:
            return _empty_handover_result()

        df = pd.DataFrame(result)

        date_column = "date_status_update"
        df[date_column] = pd.to_datetime(df[date_column], errors="coerce")
        df = df.replace([np.nan, np.inf, -np.inf], None)

        def clean_bioproject_names(values):
            unique_vals = {v for v in values if v}
            if unique_vals:
                return ", ".join(sorted(unique_vals))
            else:
                return "Unknown"

        def combine_bioproject_names(values):
            unique_vals = set()
            for value in values:
                if not value or value == "Unknown":
                    continue
                unique_vals.update(
                    part.strip() for part in str(value).split(",") if part.strip()
                )
            return ", ".join(sorted(unique_vals)) if unique_vals else "Unknown"

        def overview_grouped_status(status):
            if status in {"live", "released"}:
                return "released"
            if status in {
                "handed_over",
                "coming_soon",
            }:
                return "handed_over"
            if status in {"completed", "pre_released"}:
                return "ready_for_handover"
            if status == "in_progress":
                return "in_progress"
            if status == "abandoned":
                return "abandoned"
            if status == "archive":
                return "archive"
            if status in {
                "insufficient_data",
                "poor_genome_busco",
                "low_genome_busco",
                "check_busco",
            }:
                return "data_error"
            return "other"

        def summary_grouped_status(status):
            if status in {"live", "released"}:
                return "live"
            if status in {
                "completed",
                "complete",
                "handed_over",
                "pre_released",
                "coming_soon",
            }:
                return "pending"
            if status == "in_progress":
                return "in_progress"
            if status == "abandoned":
                return "abandoned"
            if status in {
                "insufficient_data",
                "poor_genome_busco",
                "low_genome_busco",
                "check_busco",
            }:
                return "data_error"
            return "other"

        def assign_pipeline(lowest_taxon_id):
            try:
                _, _, _, pipeline = assign_clade_and_species(
                    lowest_taxon_id, clade_data, taxonomy_dict
                )
                return pipeline
            except Exception:
                logging.warning(
                    "Could not assign pipeline for taxon_id %s", lowest_taxon_id
                )
                return "anno"

        clade_data = load_clade_data()

        collapsed = df.groupby("genebuild_status_id", as_index=False).agg(
            {
                "assembly_id": "first",
                "gca": "first",
                "bioproject_name": clean_bioproject_names,
                "date_status_update": "first",
                "gb_status": "first",
                "scientific_name": "first",
                "lowest_taxon_id": "first",
                "annotation_method": "first",
                "other_annotation_count": "max",
            }
        )
        collapsed["pipeline"] = collapsed["lowest_taxon_id"].apply(assign_pipeline)

        df_ready = collapsed[
            collapsed["gb_status"].isin(["pre_released", "completed"])
        ].copy()
        df_ready.loc[:, "core_name"] = df_ready.apply(
            lambda row: _build_core_name(row["scientific_name"], row["gca"]),
            axis=1,
        )

        # check data and in progress for more than 6 months
        date_series = pd.to_datetime(collapsed["date_status_update"], errors="coerce")
        six_months_ago = pd.Timestamp.today() - pd.DateOffset(months=6)
        collapsed["days_since_update"] = (
            pd.Timestamp.today().normalize() - date_series
        ).dt.days

        count_ho_ready = df_ready["gca"].nunique()

        # handed over (from collapsed)
        check_data = ["check_busco", "insufficient_data", "poor_genome_busco", "low_genome_busco"]
        check_data_df = collapsed[collapsed["gb_status"].isin(check_data)]

        count_data = (
            check_data_df.loc[(date_series < six_months_ago), "gca"]
        ).nunique()

        list_data = check_data_df.loc[
            (date_series < six_months_ago),
            [
                "gca",
                "scientific_name",
                "gb_status",
                "bioproject_name",
                "date_status_update",
                "days_since_update",
            ],
        ]

        # pending
        count_pending = (
            collapsed.loc[
                (collapsed["gb_status"] == "in_progress")
                & (date_series < six_months_ago),
                "gca",
            ]
        ).nunique()

        list_pending = collapsed.loc[
            (collapsed["gb_status"] == "in_progress") & (date_series < six_months_ago),
            [
                "gca",
                "scientific_name",
                "gb_status",
                "bioproject_name",
                "date_status_update",
                "days_since_update",
            ],
        ]

        overview = collapsed.copy()
        overview["dashboard_status"] = overview["gb_status"].apply(overview_grouped_status)
        overview["days_since_update"] = (
            pd.Timestamp.today().normalize() - date_series
        ).dt.days
        overview["older_than_180_days"] = overview["days_since_update"] > 180
        overview["has_other_annotation"] = overview["other_annotation_count"] > 0

        show_abandoned = (
            overview["dashboard_status"].eq("abandoned")
            & overview["older_than_180_days"]
            & ~overview["has_other_annotation"]
        )
        overview = overview.loc[
            ~overview["dashboard_status"].isin(["released", "archive"])
            & (~overview["dashboard_status"].eq("abandoned") | show_abandoned)
        ].copy()

        overview["queue"] = overview["dashboard_status"].replace(
            {
                "handed_over": "Handed over",
                "ready_for_handover": "Ready for handover",
                "abandoned": "Abandoned",
                "data_error": "Data error",
                "in_progress": "In progress",
                "other": "Other",
            }
        )
        overview["next_action"] = np.select(
            [
                overview["dashboard_status"].eq("handed_over"),
                overview["dashboard_status"].eq("ready_for_handover"),
                overview["dashboard_status"].eq("abandoned"),
                overview["dashboard_status"].eq("data_error")
                & overview["older_than_180_days"],
                overview["dashboard_status"].eq("data_error")
                & ~overview["older_than_180_days"],
                overview["dashboard_status"].eq("in_progress"),
            ],
            [
                "No action needed, waiting for release",
                "Check quality and handover",
                "Check if it can be annotated",
                "Check if it can be annotated",
                "Waiting for more data",
                "Proceed with annotation",
            ],
            default="Review current status and next action",
        )
        overview["priority"] = np.select(
            [
                overview["dashboard_status"].eq("ready_for_handover"),
                overview["dashboard_status"].eq("abandoned"),
                overview["dashboard_status"].eq("data_error")
                & overview["older_than_180_days"],
            ],
            ["high", "high", "high"],
            default="normal",
        )
        overview["days_since_update"] = overview["days_since_update"].astype(
            object
        ).where(pd.notna(overview["days_since_update"]), None)

        collapsed_for_summary = collapsed.assign(
            dashboard_status=collapsed["gb_status"].apply(summary_grouped_status),
            pipeline_normalized=collapsed["pipeline"]
            .fillna("")
            .astype(str)
            .str.lower(),
        )
        collapsed_for_summary["is_main"] = collapsed_for_summary[
            "pipeline_normalized"
        ].eq("main")
        collapsed_for_summary["is_anno"] = collapsed_for_summary[
            "pipeline_normalized"
        ].eq("anno")
        collapsed_for_summary["is_hprc"] = collapsed_for_summary[
            "pipeline_normalized"
        ].eq("hprc")
        status_summary = (
            collapsed_for_summary.groupby("dashboard_status", as_index=False)
            .agg(
                main=("is_main", "sum"),
                anno=("is_anno", "sum"),
                hprc=("is_hprc", "sum"),
                bioprojects=("bioproject_name", combine_bioproject_names),
            )
            .rename(columns={"dashboard_status": "gb_status"})
            .sort_values(["gb_status"], ascending=[True])
        )
        date_columns = ["date_status_update"]
        for column in date_columns:
            df_ready[column] = pd.to_datetime(
                df_ready[column], errors="coerce"
            ).dt.strftime("%Y-%m-%d")
            df_ready[column] = df_ready[column].where(df_ready[column].notna(), None)
            overview[column] = pd.to_datetime(
                overview[column], errors="coerce"
            ).dt.strftime("%Y-%m-%d")
            overview[column] = overview[column].where(overview[column].notna(), None)

        def stale_records(frame):
            if frame.empty:
                return []
            frame = frame.copy()
            frame["date_status_update"] = pd.to_datetime(
                frame["date_status_update"], errors="coerce"
            ).dt.strftime("%Y-%m-%d")
            frame["date_status_update"] = frame["date_status_update"].where(
                frame["date_status_update"].notna(), None
            )
            frame["days_since_update"] = frame["days_since_update"].astype(
                object
            ).where(pd.notna(frame["days_since_update"]), None)
            return _json_safe_records(frame)

        return (
            _json_safe_records(
                df_ready[
                    [
                        "gca",
                        "scientific_name",
                        "bioproject_name",
                        "date_status_update",
                        "core_name",
                    ]
                ]
            ),
            int(count_ho_ready),
            int(count_data),
            int(count_pending),
            stale_records(list_data),
            stale_records(list_pending),
            _json_safe_records(overview[
                [
                    "gca",
                    "genebuild_status_id",
                    "scientific_name",
                    "gb_status",
                    "dashboard_status",
                    "queue",
                    "next_action",
                    "priority",
                    "bioproject_name",
                    "annotation_method",
                    "date_status_update",
                    "days_since_update",
                ]
            ]),
            _json_safe_records(status_summary),
        )

    except Exception as e:
        logging.error(f"Error in get_ready_to_ho: {e}", exc_info=True)
        return {}


def update_gca(genebuilder, items, new_status):
    try:
        if not items:
            return 0

        with get_db_connection("meta_write") as conn:
            cursor = conn.cursor(pymysql.cursors.DictCursor)

            query = """
				UPDATE
					genebuild_status g
				SET g.gb_status = %s,
                    g.date_status_update = CURRENT_DATE
				WHERE g.genebuilder = %s
				  AND g.gca_accession = %s
				  AND g.annotation_method = %s
			"""
            updated_rows = 0

            for item in items:
                cursor.execute(
                    query,
                    (
                        new_status,
                        genebuilder,
                        item["gca"],
                        item["annotation_method"],
                    ),
                )
                updated_rows += cursor.rowcount
                if cursor.rowcount == 0:
                    logging.warning(
                        f"No rows updated for GCA={item['gca']} "
                        f"method={item['annotation_method']} "
                        f"genebuilder={genebuilder}"
                        f"new status={new_status}"
                    )

            conn.commit()
            return updated_rows

    except Exception as e:
        logging.error(f"Error in update_gca: {e}", exc_info=True)
        return 0
