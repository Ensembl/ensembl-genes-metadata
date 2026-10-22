# app/services/home_page_service.py
import datetime
import logging

import pandas as pd

from metadata_app.backend.app.core.database import get_db_connection
from metadata_app.backend.app.services.filter_scope_service import add_scope_conditions


def get_annotation_counts_by_bioproject():
    """Returns a count of annotations per BioProject in main_bioproject, with annotation and assembly info."""
    try:
        with get_db_connection("meta") as conn:
            cursor = conn.cursor()
            query = """
                SELECT 
                    b.bioproject_id,
                    mb.bioproject_name,
                    COUNT(g.assembly_id) AS annotation_count,

                    (
                        SELECT COUNT(DISTINCT a.assembly_id)
                        FROM assembly a
                        LEFT JOIN bioproject b2 ON a.assembly_id = b2.assembly_id
                        LEFT JOIN genebuild_status g2 ON a.assembly_id = g2.assembly_id

                        WHERE b2.bioproject_id = b.bioproject_id
                          AND g2.assembly_id IS NULL
                          AND a.asm_name NOT LIKE "%alternate_haplotype%"
                          AND a.asm_level IN ('Chromosome', 'Complete genome')
                          AND a.is_current = 'current'
                    ) AS qualified_assembly_count,

                    (
                        SELECT COUNT(DISTINCT g2.assembly_id)
                        FROM genebuild_status g2
                        JOIN bioproject b3 ON g2.assembly_id = b3.assembly_id
                        WHERE g2.gb_status IN ('in_progress', 'complete', 'pre_released', 'coming_soon')
                          AND b3.bioproject_id = b.bioproject_id
                    ) AS in_progress

                FROM bioproject b
                JOIN genebuild_status g ON b.assembly_id = g.assembly_id
                JOIN main_bioproject mb ON mb.bioproject_id = b.bioproject_id
                WHERE g.gb_status = 'live'
                GROUP BY b.bioproject_id, mb.bioproject_name
                
            """
            cursor.execute(query)
            result = cursor.fetchall()

        df = pd.DataFrame(
            result,
            columns=[
                "bioproject_id",
                "bioproject_name",
                "annotation_count",
                "qualified_assembly_count",
                "in_progress",
            ],
        )
        return df.to_dict(orient="records")

    except Exception as e:
        logging.error(f"Error fetching annotation counts: {e}")
        return []


def get_annotation_counts_by_group():
    """Returns a count of annotations per BioProject in main_bioproject, with annotation and assembly info."""
    try:
        with get_db_connection("meta") as conn:
            cursor = conn.cursor()
            query = """
                   SELECT
					    cg.group_id,
					    cg.group_name,
					    COUNT(DISTINCT CASE WHEN g.gb_status = 'live' THEN a.assembly_id END) AS annotation_count,
					    COUNT(DISTINCT CASE
					        WHEN g.assembly_id IS NULL
					             AND a.asm_name NOT LIKE "%alternate_haplotype%"
					             AND a.asm_level IN ('Chromosome', 'Complete genome')
					             AND a.is_current = 'current'
					        THEN a.assembly_id
					    END) AS qualified_assembly_count,
					    COUNT(DISTINCT CASE
					        WHEN g.gb_status IN ('in_progress', 'complete', 'pre_released', 'coming_soon')
					        THEN a.assembly_id
					    END) AS in_progress
					FROM assembly a
					JOIN custom_group cg
					    ON (cg.group_type = 'taxon' AND a.lowest_taxon_id = cg.item)
					    OR (cg.group_type = 'assembly' AND a.gca_chain = cg.item)
					LEFT JOIN genebuild_status g ON a.assembly_id = g.assembly_id
					GROUP BY cg.group_name;
			        """
            cursor.execute(query)
            result = cursor.fetchall()

        df = pd.DataFrame(
            result,
            columns=[
                "group_id",
                "group_name",
                "annotation_count",
                "qualified_assembly_count",
                "in_progress",
            ],
        )
        return df.to_dict(orient="records")

    except Exception as e:
        logging.error(f"Error fetching annotation counts: {e}")
        return []


def get_assemblies_per_year():
    try:
        with get_db_connection("meta") as conn:
            cursor = conn.cursor()
            query = """
                SELECT a.assembly_id, a.release_date
                FROM assembly a
                WHERE a.is_current = 'current'
            """
            cursor.execute(query)
            result = cursor.fetchall()

        df = pd.DataFrame(result, columns=["assembly_id", "release_date"])
        if df.empty:
            return []

        # Convert release_date to datetime and extract year
        df["release_date"] = pd.to_datetime(df["release_date"], errors="coerce")
        df["year"] = df["release_date"].dt.year

        # Filter to include only 2019 and later
        df = df[df["year"] >= 2019]

        # Group by year, then carry the annual totals forward cumulatively.
        counts_by_year = df.groupby("year")["assembly_id"].nunique().reset_index().sort_values("year")
        counts_by_year.columns = ["year", "assembly_count"]
        counts_by_year["assembly_count"] = counts_by_year["assembly_count"].cumsum()

        return counts_by_year.to_dict(orient="records")

    except Exception as e:
        print(f"Error in get_assemblies_per_year: {e}")
        return []


def get_annotations_per_year():
    try:
        with get_db_connection("meta") as conn:
            cursor = conn.cursor()
            query = """
                SELECT g.genebuild_status_id, g.date_status_update
                FROM genebuild_status g
            """
            cursor.execute(query)
            result = cursor.fetchall()

        df = pd.DataFrame(result, columns=["genebuild_status_id", "date_status_update"])
        if df.empty:
            return []

        # Convert date_completed to datetime and extract year
        df["date_status_update"] = pd.to_datetime(df["date_status_update"], errors="coerce")
        df["year"] = df["date_status_update"].dt.year
        # Filter to include only 2019 and later
        df = df[df["year"] >= 2019]

        # Group by year and count
        counts_by_year = df.groupby("year")["genebuild_status_id"].nunique().reset_index().sort_values("year")
        counts_by_year.columns = ["year", "annotation_count"]
        counts_by_year["annotation_count"] = counts_by_year["annotation_count"].cumsum()

        return counts_by_year.to_dict(orient="records")

    except Exception as e:
        print(f"Error in get_annotations_per_year: {e}")
        return []


def get_annotation_growth_by_quarter():
    """Return distinct annotated genomes by quarter for the last four years."""
    try:
        with get_db_connection("meta") as conn:
            cursor = conn.cursor()
            query = """
                SELECT
                    YEAR(g.release_date) AS year,
                    QUARTER(g.release_date) AS quarter_number,
                    g.annotation_method,
                    COUNT(DISTINCT g.genebuild_status_id) AS genome_count
                FROM genebuild_status g
                WHERE g.release_date >= DATE_SUB(CURDATE(), INTERVAL 4 YEAR)
                  AND g.annotation_method IN ('full_genebuild', 'anno', 'braker', 'helixer')
                  AND g.gb_status = 'live'
                GROUP BY year, quarter_number, annotation_method
                ORDER BY year, quarter_number, annotation_method
            """
            cursor.execute(query)
            result = cursor.fetchall()

            baseline_query = """
                SELECT g.annotation_method, COUNT(DISTINCT g.genebuild_status_id) AS genome_count
                FROM genebuild_status g
                WHERE g.release_date < DATE_SUB(CURDATE(), INTERVAL 4 YEAR)
                  AND g.annotation_method IN ('full_genebuild', 'anno', 'braker', 'helixer')
                  AND g.gb_status = 'live'
                GROUP BY g.annotation_method
            """
            cursor.execute(baseline_query)
            baseline_result = cursor.fetchall()

        columns = ["year", "quarter_number", "annotation_method", "genome_count"]
        df = pd.DataFrame(result, columns=columns)
        if df.empty:
            wide = pd.DataFrame(columns=["year", "quarter_number", "quarter"])
        else:
            df["quarter"] = df["year"].astype(str) + " Q" + df["quarter_number"].astype(str)
            wide = (
                df.pivot_table(
                    index=["year", "quarter_number", "quarter"],
                    columns="annotation_method",
                    values="genome_count",
                    aggfunc="sum",
                    fill_value=0,
                )
                .reset_index()
                .sort_values(["year", "quarter_number"])
            )
        start_period = (pd.Timestamp.today() - pd.DateOffset(years=4)).to_period("Q")
        end_period = pd.Timestamp.today().to_period("Q")
        periods = pd.period_range(start=start_period, end=end_period, freq="Q")
        quarter_index = pd.DataFrame(
            {
                "year": periods.year,
                "quarter_number": periods.quarter,
                "quarter": [f"{period.year} Q{period.quarter}" for period in periods],
            }
        )
        wide = quarter_index.merge(wide, on=["year", "quarter_number", "quarter"], how="left")
        methods = ["full_genebuild", "anno", "braker", "helixer"]
        for method in methods:
            if method not in wide:
                wide[method] = 0
            wide[method] = wide[method].fillna(0).astype(int)
        wide[methods] = wide[methods].cumsum()
        baseline = dict(baseline_result)
        for method in methods:
            wide[method] += int(baseline.get(method, 0))
        wide["total_genomes"] = wide[methods].sum(axis=1)
        return wide[["quarter", "total_genomes", *methods]].to_dict(orient="records")

    except Exception as e:
        logging.error(f"Error fetching quarterly annotation growth: {e}")
        return []


def get_metadata_registry_update_dates():
    try:
        with get_db_connection("meta") as conn:
            cursor = conn.cursor()
            query = """
                SELECT date_value
                FROM update_date
                WHERE update_type = 'regular_update'
            """
            cursor.execute(query)
            result = cursor.fetchall()

        # Convert result to DataFrame with correct column name
        df = pd.DataFrame(result, columns=["date_value"])

        if df.empty:
            return []

        # Convert to string list of dates
        return df["date_value"].astype(str).tolist()

    except Exception as e:
        print(f"Error fetching metadata registry update dates: {e}")
        return []


def get_transcriptomic_registry_update_dates():
    try:
        with get_db_connection("transcriptomic") as conn:
            cursor = conn.cursor()
            query = """
                SELECT last_check
                FROM meta
            """
            cursor.execute(query)
            result = cursor.fetchall()

        # Convert result to DataFrame with correct column name
        df = pd.DataFrame(result, columns=["last_check"])

        if df.empty:
            return []

        # Convert to string list of dates
        return df["last_check"].astype(str).tolist()

    except Exception as e:
        print(f"Error fetching transcriptomic registry update dates: {e}")
        return []


def bin_by_genebuild_method(bioproject_id, taxon_id, release_date):
    """Bins assemblies based on genebuild.method."""
    try:
        with get_db_connection("meta") as conn:
            cursor = conn.cursor()
            # Build query conditions
            conditions = []
            params = []

            add_scope_conditions(
                conditions,
                params,
                bioproject_id=bioproject_id,
                taxon_id=taxon_id,
                taxon_column="a.lowest_taxon_id",
            )

            if release_date:
                if isinstance(release_date, pd.Timestamp):
                    release_date = release_date.strftime("%Y-%m-%d")
                elif isinstance(release_date, (datetime.date, datetime.datetime)):
                    release_date = release_date.strftime("%Y-%m-%d")
                conditions.append("g.release_date_beta >= %s")
                params.append(release_date)
                logging.info(f"Filtering by release date: {release_date}")

            # If there are conditions, join them with AND; otherwise, select all
            where_clause = " WHERE " + " AND ".join(conditions) if conditions else ""

            query = f"""
                    SELECT g.annotation_method, g.release_date, b.bioproject_id, a.lowest_taxon_id, g.genebuild_status_id
                    FROM genebuild_status g
                    JOIN assembly a ON g.assembly_id = a.assembly_id
                    JOIN bioproject b on g.assembly_id = b.assembly_id
                    {where_clause};
                """
            cursor.execute(query, params)
            result = cursor.fetchall()

        df = pd.DataFrame(
            result,
            columns=[
                "genebuild_status_id",
                "annotation_method",
                "release_date",
                "bioproject_id",
                "lowest_taxon_id",
            ],
        )

        df = df.drop_duplicates(subset="genebuild_status_id", keep="first")

        if df.empty:
            return []

        # Group by genebuilder and count
        method_summary = (
            df.groupby("annotation_method", observed=False).size().reset_index(name="number_of_annotations")
        )

        return method_summary.to_dict(orient="records")

    except Exception as e:
        # Optionally log or re-raise
        print(f"Error in bin_by_genebuld_method: {e}")
        return []


def annotation_time():
    """Return recent duration summaries and the historical monthly trend."""
    try:
        with get_db_connection("meta") as conn:
            cursor = conn.cursor()
            query = """
            SELECT
                g.genebuild_status_id,
                g.annotation_method,
                g.date_started,
                g.last_genebuild_update AS last_geneset_update
            FROM genebuild_status g
            WHERE g.annotation_method IN ('anno', 'full_genebuild')
              AND g.gb_status = 'live'  
              AND g.date_started IS NOT NULL
            """
            cursor.execute(query)
            result = cursor.fetchall()

        columns = ["genebuild_status_id", "annotation_method", "date_started", "last_geneset_update"]
        df = pd.DataFrame(result, columns=columns)
        if df.empty:
            return {"summary": [], "trend": [], "window_months": 3}

        df["date_started"] = pd.to_datetime(df["date_started"], errors="coerce")
        df["last_geneset_update"] = pd.to_datetime(df["last_geneset_update"], errors="coerce")
        df["duration_days"] = (df["last_geneset_update"] - df["date_started"]).dt.total_seconds() / (
            24 * 60 * 60
        )
        df = df.dropna(subset=["annotation_method", "last_geneset_update", "duration_days"])
        df = df[(df["duration_days"] >= 0) & (df["duration_days"] <= 3650)]
        df = df.drop_duplicates(subset=["genebuild_status_id", "annotation_method"], keep="first")
        if df.empty:
            return {"summary": [], "trend": [], "window_months": 3}

        window_end = df["last_geneset_update"].max()
        window_start = window_end - pd.DateOffset(months=3)
        recent = df[df["last_geneset_update"] >= window_start]
        methods = ["anno", "full_genebuild", "projection_build"]
        summary = []
        for method in methods:
            durations = recent.loc[recent["annotation_method"] == method, "duration_days"]
            summary.append(
                {
                    "annotation_method": method,
                    "sample_size": int(durations.count()),
                    "average_days": round(float(durations.mean()), 1) if not durations.empty else None,
                    "median_days": round(float(durations.median()), 1) if not durations.empty else None,
                    "minimum_days": round(float(durations.min()), 1) if not durations.empty else None,
                    "maximum_days": round(float(durations.max()), 1) if not durations.empty else None,
                }
            )

        df["month"] = df["last_geneset_update"].dt.to_period("M").astype(str)
        trend = (
            df.groupby(["month", "annotation_method"], as_index=False)
            .agg(average_days=("duration_days", "mean"), sample_size=("duration_days", "count"))
            .pivot(index="month", columns="annotation_method", values=["average_days", "sample_size"])
            .reset_index()
        )
        trend.columns = [
            "month" if column[0] == "month" else f"{column[0]}_{column[1]}"
            for column in trend.columns.to_flat_index()
        ]
        for method in methods:
            for prefix in ("average_days", "sample_size"):
                column = f"{prefix}_{method}"
                if column not in trend:
                    trend[column] = None
        trend = trend.sort_values("month")
        for column in trend.columns:
            if column.startswith("average_days_"):
                trend[column] = pd.to_numeric(trend[column], errors="coerce").round(1)
            elif column.startswith("sample_size_"):
                trend[column] = trend[column].fillna(0).astype(int)
        trend = trend.astype(object).where(pd.notna(trend), None)

        return {
            "summary": summary,
            "trend": trend.to_dict(orient="records"),
            "window_months": 3,
            "window_start": window_start.strftime("%Y-%m"),
            "window_end": window_end.strftime("%Y-%m"),
        }

    except Exception as e:
        logging.error(f"Error fetching annotation durations: {e}", exc_info=True)
        return {"summary": [], "trend": [], "window_months": 3}
