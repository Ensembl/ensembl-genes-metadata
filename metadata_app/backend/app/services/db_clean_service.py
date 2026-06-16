# app/services/db_clean_service.py

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
from metadata_app.backend.app.core.database import get_db_connection
import pymysql.cursors
import pymysql
import pandas as pd
import re

# Example server info: list of dicts
servers = [
    {
        "host": "mysql-ens-genebuild-prod-1",
        "port": 4527,
        "user": "ensro",
        "password": "",
    },
    {
        "host": "mysql-ens-genebuild-prod-2",
        "port": 4528,
        "user": "ensro",
        "password": "",
    },
    {
        "host": "mysql-ens-genebuild-prod-3",
        "port": 4529,
        "user": "ensro",
        "password": "",
    },
    {
        "host": "mysql-ens-genebuild-prod-4",
        "port": 4530,
        "user": "ensro",
        "password": "",
    },
    {
        "host": "mysql-ens-genebuild-prod-5",
        "port": 4531,
        "user": "ensro",
        "password": "",
    },
    {
        "host": "mysql-ens-genebuild-prod-6",
        "port": 4532,
        "user": "ensro",
        "password": "",
    },
    {
        "host": "mysql-ens-genebuild-prod-7",
        "port": 4533,
        "user": "ensro",
        "password": "",
    },
]


def get_live_annotations(genebuilder):
    try:
        with get_db_connection("meta") as conn:
            with conn.cursor(pymysql.cursors.DictCursor) as cursor:
                query = """
                    SELECT 
                        g.genebuild_status_id,
                        g.gb_status,
                        g.last_genebuild_update,
                        g.genebuilder,
                        s.scientific_name,
                        CONCAT(a.gca_chain, ".", a.gca_version) AS gca
                    FROM genebuild_status g
                    JOIN assembly a ON a.assembly_id = g.assembly_id
                    LEFT JOIN species s ON s.lowest_taxon_id = a.lowest_taxon_id
                    WHERE g.gb_status = 'live'
                    AND g.genebuilder = %s
                """
                cursor.execute(query, (genebuilder,))
                result = cursor.fetchall()

        # Convert to DataFrame
        df = pd.DataFrame(result)

        # Replace NaN/Inf/-Inf with "Unknown"
        df = df.replace([np.nan, np.inf, -np.inf], "Unknown")

        return df

    except Exception as e:
        logging.error(
            f"Error fetching live annotations for {genebuilder}: {e}", exc_info=True
        )
        return pd.DataFrame()


def _normalise_gca_for_db_name(gca):
    return str(gca).lower().replace("gca_", "gca").replace(".", "v")


def _scientific_name_to_db_prefix(scientific_name):
    if pd.isna(scientific_name):
        return None

    prefix = re.sub(r"\s+", "_", str(scientific_name).strip().lower())
    return prefix or None


def _build_database_matchers(df, chunk_size=250):
    """
    Compile database-name matchers once instead of rebuilding/checking one
    regex per live annotation for every database on every server.
    """
    candidates = set()

    for _, row in df.iterrows():
        gca_clean = _normalise_gca_for_db_name(row["gca"])
        prefixes = {str(row["genebuilder"])}

        if "scientific_name" in row:
            species_prefix = _scientific_name_to_db_prefix(row["scientific_name"])
            if species_prefix:
                prefixes.add(species_prefix)

        for prefix in prefixes:
            candidates.add(f"{prefix}_{gca_clean}")

    escaped_candidates = [re.escape(candidate) for candidate in candidates]
    return [
        re.compile(rf"^(?:{'|'.join(escaped_candidates[i:i + chunk_size])}).*$", re.IGNORECASE)
        for i in range(0, len(escaped_candidates), chunk_size)
    ]


def _database_matches(db_name, matchers):
    return any(matcher.fullmatch(db_name) for matcher in matchers)


def _find_databases_on_server(server, matchers, genebuilder):
    db_list = []
    conn = None

    try:
        conn = pymysql.connect(
            host=server["host"],
            port=server["port"],
            user=server["user"],
            password=server["password"],
            cursorclass=pymysql.cursors.DictCursor,
            connect_timeout=5,
            read_timeout=30,
        )
        with conn.cursor() as cursor:
            cursor.execute("SHOW DATABASES;")
            databases = cursor.fetchall()
            for db in databases:
                db_name = db["Database"]
                if _database_matches(db_name, matchers):
                    db_list.append(
                        {
                            "database": db_name,
                            "server": server["host"],
                            "port": server["port"],
                            "genebuilder": genebuilder,
                        }
                    )
    except Exception as e:
        logging.error(
            f"Error connecting to {server['host']}:{server['port']} - {e}",
            exc_info=True,
        )
    finally:
        if conn:
            conn.close()

    return db_list


def find_genebuilder_databases(df):
    """
    Find genebuilder-related databases across servers based on live annotations.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain columns "gca" and "genebuilder".
        If "scientific_name" is present, species-prefixed database names are
        also matched for GB1-style core database names.

    Returns
    -------
    pd.DataFrame
        Columns: database, server, port, genebuilder
    """
    if df.empty:
        return pd.DataFrame(columns=["database", "server", "port"])

    db_list = []
    genebuilder = str(df.iloc[0]["genebuilder"])
    matchers = _build_database_matchers(df)

    with ThreadPoolExecutor(max_workers=len(servers)) as executor:
        futures = [
            executor.submit(_find_databases_on_server, server, matchers, genebuilder)
            for server in servers
        ]

        for future in as_completed(futures):
            db_list.extend(future.result())

    if not db_list:
        return pd.DataFrame(columns=["database", "server", "port"])

    df_result = pd.DataFrame(db_list)

    return df_result


def generate_drop_script(df):
    """
    Generate a MySQL-executable script to drop databases, grouped by server/port.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain 'database', 'server', 'port'

    Returns
    -------
    str
        SQL script suitable for direct execution via mysql client.
    """
    if df.empty:
        return "-- No databases found to drop.\n"

    script_lines = [
        "-- Cleanup script for live databases",
        "-- ========================================",
        "-- Generated automatically by db_clean_service",
        "-- Please be careful when using this script. Anno pipelines cannot be checked per live GCA. Check if they can be deleted.",
        "",
    ]

    # Group by server/port to provide connection context
    for (server, port), group in df.groupby(["server", "port"]):
        script_lines.append(f"-- Server: {server}:{port}")
        script_lines.append(f"-- To execute on this server:")
        script_lines.append(f"-- mysql -h{server} -P{port} -uensadmin -p")
        script_lines.append("")  # blank line
        for _, row in group.iterrows():
            db_name = row["database"]
            script_lines.append(f"DROP DATABASE IF EXISTS `{db_name}`;")
        script_lines.append("")  # blank line between servers

    return "\n".join(script_lines)


def server_clean_main(genebuilder):
    df = get_live_annotations(genebuilder)
    df_db = find_genebuilder_databases(df)
    commands = generate_drop_script(df_db)

    df_db = df_db.to_dict(orient="records")
    return df_db, commands
