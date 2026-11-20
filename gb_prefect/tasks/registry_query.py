"""Collect a GCA list from the metadata registry with a SQL query."""

import re
from pathlib import Path
from typing import Any, Dict, List

import pymysql  # type: ignore
from prefect import task  # type: ignore

# The default list: every GCA in genebuild_status whose genebuild isn't finished.
DEFAULT_GCA_QUERY = (
    "SELECT gca_accession FROM genebuild_status "
    "WHERE gb_status NOT IN ('handed_over', 'live', 'coming_soon')"
)

GCA_PATTERN = re.compile(r"^GC[AF]_\d{9}\.\d+$")


def check_select_query(query: str) -> str:
    """Return query without surrounding whitespace or a trailing ';' after checking it is one
    SELECT (or WITH ... SELECT) statement. The query also runs in a read-only transaction,
    so a check missed here still can't change the DB."""
    cleaned = query.strip().rstrip(";").strip()
    if ";" in cleaned:
        raise ValueError(f"Give a single SQL statement, not several: {query}")
    if not re.match(r"^(SELECT|WITH)\b", cleaned, re.IGNORECASE):
        raise ValueError(f"Only SELECT queries are allowed: {query}")
    return cleaned


def unique_gcas(rows: List[Any]) -> List[str]:
    """Return the first column of rows as accessions, de-duplicated in query order.

    Raises ValueError if a value doesn't look like an accession (e.g. the query selects the
    wrong column first)."""
    gcas: List[str] = []
    seen = set()
    for row in rows:
        value = str(row[0]).strip() if row and row[0] is not None else ""
        if not GCA_PATTERN.match(value):
            raise ValueError(f"The query's first column must be a GCA accession, got: {value!r}")
        if value not in seen:
            seen.add(value)
            gcas.append(value)
    return gcas


def run_select_query(query: str, db_params: Dict[str, Any]) -> List[Any]:
    """Run query in a read-only transaction and return every row."""
    with pymysql.connect(**db_params) as conn:
        with conn.cursor() as cursor:
            cursor.execute("START TRANSACTION READ ONLY")
            cursor.execute(query)
            rows = list(cursor.fetchall())
        conn.rollback()
    return rows


def write_gca_file(output_file: str, gcas: List[str]) -> str:
    """Write gcas to output_file, one per line (creating its folder). Returns output_file."""
    path = Path(output_file)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"{gca}\n" for gca in gcas), encoding="utf-8")
    return output_file


def split_batches(gcas: List[str], batch_size: int) -> List[List[str]]:
    """Split gcas, in order, into batches of batch_size (the last one may be smaller)."""
    if batch_size < 1:
        raise ValueError(f"batch_size must be at least 1, got {batch_size}")
    return [gcas[start : start + batch_size] for start in range(0, len(gcas), batch_size)]


@task(log_prints=True)
def collect_gca_list(query: str, db_params: Dict[str, Any], output_file: str) -> List[str]:
    """Run a SELECT query on the metadata DB and write the GCAs it returns (its first column,
    de-duplicated) to output_file, one per line. Returns the GCAs."""
    query = check_select_query(query)
    print(f"Collecting GCAs with:\n{query}")
    gcas = unique_gcas(run_select_query(query, db_params))
    write_gca_file(output_file, gcas)
    print(f"Wrote {len(gcas)} GCA(s) to {output_file}")
    return gcas
