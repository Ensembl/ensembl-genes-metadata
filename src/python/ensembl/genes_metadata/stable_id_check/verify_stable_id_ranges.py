#!/usr/bin/env python3
"""Verify stable ID ranges in core databases against the metadata registry."""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from pathlib import Path
from typing import Any


def load_driver():
    try:
        import pymysql

        return pymysql
    except ImportError:
        pass
    try:
        import mysql.connector

        return mysql.connector
    except ImportError as exc:
        raise RuntimeError("install pymysql or mysql-connector-python") from exc


def connect_server(database: str) -> Any:
    driver = load_driver()
    kwargs: dict[str, Any] = {
        "host": "mysql-ens-genebuild-prod-1",
        "port": 4527,
        "user": "ensro",
        "database": database,
    }

    return driver.connect(**kwargs)


def connect_genebuild(core_db: str) -> Any:
    driver = load_driver()
    kwargs: dict[str, Any] = {
        "host": "mysql-ens-genebuild-prod-1",
        "port": 4527,
        "user": "ensro",
        "database": core_db,
    }

    return driver.connect(**kwargs)


def connect_metadata_registry() -> Any:
    driver = load_driver()
    kwargs: dict[str, Any] = {
        "host": "mysql-ens-genebuild-prod-1",
        "port": 4527,
        "user": "ensro",
        "database": "gb_assembly_metadata",
    }

    return driver.connect(**kwargs)


def parse_stable_id(stable_id: str) -> tuple[str, int]:
    match = re.match(r"^([A-Za-z]+)(\d+)$", stable_id)
    if not match:
        raise ValueError(f"cannot parse stable_id: {stable_id}")
    return match.group(1), int(match.group(2))


def query_core(core_db: str) -> dict[str, Any]:
    conn = connect_genebuild(core_db)
    try:
        cursor = conn.cursor()
        cursor.execute("SELECT meta_value FROM meta WHERE meta_key = 'assembly.accession'")
        row = cursor.fetchone()
        gca_accession = row[0] if row else None

        cursor.execute("SELECT MIN(stable_id), MAX(stable_id) FROM gene")
        min_id, max_id = cursor.fetchone()

        prefix_min, num_min = parse_stable_id(min_id)
        prefix_max, num_max = parse_stable_id(max_id)

        if prefix_min != prefix_max:
            raise ValueError(f"prefix mismatch: min={prefix_min}, max={prefix_max}")

        return {
            "gca_accession": gca_accession,
            "prefix": prefix_min,
            "min_stable_id": min_id,
            "max_stable_id": max_id,
            "min_number": num_min,
            "max_number": num_max,
        }
    finally:
        cursor.close()
        conn.close()


def query_registry(gca_accession: str) -> dict[str, Any]:
    conn = connect_metadata_registry()
    try:
        cursor = conn.cursor()

        gca_chain, gca_version = gca_accession.rsplit(".", 1)
        cursor.execute(
            """
            SELECT * FROM species_prefix
            JOIN assembly a ON species_prefix.lowest_taxon_id = a.lowest_taxon_id
            WHERE gca_chain = %s AND gca_version = %s
            """,
            (gca_chain, int(gca_version)),
        )
        prefix_rows = cursor.fetchall()
        if prefix_rows:
            columns = [col[0] for col in cursor.description]
            prefix_rows = [dict(zip(columns, row)) for row in prefix_rows]

        cursor.execute(
            """
            SELECT * FROM stable_space
            JOIN gb_assembly_metadata.stable_space_species_log sssl
                ON stable_space.stable_space_id = sssl.stable_space_id
            WHERE gca_accession = %s
            """,
            (gca_accession,),
        )
        space_row = cursor.fetchone()
        if space_row:
            columns = [col[0] for col in cursor.description]
            space_row = dict(zip(columns, space_row))

        return {"prefix_rows": prefix_rows, "space_row": space_row}
    finally:
        cursor.close()
        conn.close()


def compare(core_info: dict[str, Any], registry_info: dict[str, Any]) -> list[str]:
    issues: list[str] = []
    prefix_rows = registry_info.get("prefix_rows", [])
    space_row = registry_info.get("space_row")

    if not prefix_rows:
        issues.append("no prefix row found in species_prefix for this GCA")
    elif len(prefix_rows) > 1:
        issues.append(f"multiple prefix rows in registry ({len(prefix_rows)})")
        for i, row in enumerate(prefix_rows, 1):
            issues.append(
                f"  registry row {i}: prefix={row.get('prefix')}, "
                f"lowest_taxon_id={row.get('lowest_taxon_id')}"
            )
        issues.append(
            f"  core values: prefix={core_info['prefix']}, "
            f"min={core_info['min_stable_id']}, max={core_info['max_stable_id']}"
        )
    else:
        reg_prefix = prefix_rows[0].get("prefix")
        core_prefix = core_info["prefix"]
        if core_prefix.endswith("G"):
            core_prefix = core_prefix[:-1]
        if reg_prefix and reg_prefix != core_prefix:
            issues.append(f"prefix issue: core={core_info['prefix']} (stripped={core_prefix}), registry={reg_prefix}")

    if not space_row:
        issues.append("no stable_space row found for this GCA")
    else:
        reg_start = space_row.get("stable_space_start")
        reg_end = space_row.get("stable_space_end")
        core_min = core_info["min_number"]
        core_max = core_info["max_number"]

        if reg_start is not None and core_min < reg_start:
            issues.append(f"range issue: core min {core_min} < registry start {reg_start}")
        if reg_end is not None and core_max > reg_end:
            issues.append(f"range issue: core max {core_max} > registry end {reg_end}")

    return issues


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("core_dbs", nargs="*", help="core database names to check")
    parser.add_argument("--db-list", type=str, help="file with one core database name per line")
    parser.add_argument("--out", type=Path, help="write results to this CSV file")
    return parser


def main(argv: list[str]) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    core_dbs: list[str] = list(args.core_dbs) if args.core_dbs else []
    if args.db_list:
        with open(args.db_list) as f:
            core_dbs.extend(line.strip() for line in f if line.strip() and not line.startswith("#"))

    if not core_dbs:
        raise SystemExit("provide core database names as arguments or via --db-list")

    print(f"{'Core DB':<40} {'GCA':<25} {'Prefix':<15} {'Min ID':<30} {'Max ID':<30} {'Status'}")
    print("-" * 160)

    rows: list[dict[str, str]] = []
    for core_db in core_dbs:
        try:
            core_info = query_core(core_db)
            registry_info = query_registry(core_info["gca_accession"])
            issues = compare(core_info, registry_info)

            status = "OK" if not issues else "; ".join(issues)
            print(
                f"{core_db:<40} {core_info['gca_accession'] or 'N/A':<25} "
                f"{core_info['prefix']:<15} {core_info['min_stable_id']:<30} "
                f"{core_info['max_stable_id']:<30} {status}"
            )
            rows.append({
                "core_db": core_db,
                "gca_accession": core_info.get("gca_accession") or "",
                "prefix": core_info["prefix"],
                "min_stable_id": core_info["min_stable_id"],
                "max_stable_id": core_info["max_stable_id"],
                "min_number": str(core_info["min_number"]),
                "max_number": str(core_info["max_number"]),
                "status": status,
            })
        except Exception as exc:
            print(f"{core_db:<40} ERROR: {exc}")
            rows.append({
                "core_db": core_db,
                "gca_accession": "",
                "prefix": "",
                "min_stable_id": "",
                "max_stable_id": "",
                "min_number": "",
                "max_number": "",
                "status": f"ERROR: {exc}",
            })

    if args.out:
        fieldnames = ["core_db", "gca_accession", "prefix", "min_stable_id", "max_stable_id", "min_number", "max_number", "status"]
        with args.out.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nCSV written to {args.out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
