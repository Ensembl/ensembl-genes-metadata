#!/usr/bin/env python3

#  See the NOTICE file distributed with this work for additional information
#  regarding copyright ownership.
#
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#  http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

"""
Given a list of GCA accessions, query the NCBI Datasets API to determine
whether each assembly is labeled as the reference assembly for its taxon.

Usage
-----
    # Pass accessions directly on the command line:
    python check_reference.py GCA_004027535.1 GCA_000001405.29

    # Or supply a file (one GCA per line, or a single-column CSV):
    python check_reference.py --file gcas.txt
    python check_reference.py --file gcas.csv

    # Choose a custom output path (default: reference_check.csv):
    python check_reference.py --file gcas.txt --output my_results.csv

Output columns
--------------
    gca               - input accession
    species_name      - organism_name from NCBI
    taxon_id          - NCBI taxonomy ID
    is_reference      - True / False / No reference available
    current_reference - GCA (or GCF paired) accession of the reference,
                        or "No reference available"
"""

import argparse
import csv
import sys
import time
from pathlib import Path

import requests

NCBI_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"
RETRY_WAIT = 2   # seconds between retries
MAX_RETRIES = 3


# ---------------------------------------------------------------------------
# API helpers
# ---------------------------------------------------------------------------

def _get(url: str) -> dict:
    """GET a URL with simple retry logic; raise on persistent failure."""
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            resp = requests.get(url, timeout=30)
            resp.raise_for_status()
            return resp.json()
        except requests.RequestException as exc:
            if attempt == MAX_RETRIES:
                raise
            print(f"  [warn] attempt {attempt} failed ({exc}); retrying …", file=sys.stderr)
            time.sleep(RETRY_WAIT)


def get_assembly_info(accession: str) -> dict:
    """
    Return a dict with keys: taxon_id, species_name.
    Raises if the accession is not found.
    """
    url = f"{NCBI_BASE}/genome/accession/{accession}/dataset_report?filters.exclude_atypical=false&filters.assembly_version=all_assemblies"
    data = _get(url)

    reports = data.get("reports", [])
    if not reports:
        raise ValueError(f"No report found for accession {accession!r}")

    organism = reports[0].get("organism", {})
    return {
        "taxon_id": organism.get("tax_id", ""),
        "species_name": organism.get("organism_name", ""),
    }


def get_reference_for_taxon(taxon_id: str) -> dict | None:
    """
    Return the first reference-assembly report for taxon_id, or None if
    no reference exists.
    """
    url = (
        f"{NCBI_BASE}/genome/taxon/{taxon_id}/dataset_report"
        f"?filters.reference_only=true"
    )
    data = _get(url)

    total = data.get("total_count", 0)
    reports = data.get("reports", [])
    if total == 0 or not reports:
        return None
    return reports[0]


def check_assembly(accession: str) -> dict:
    """
    Return a result row dict for one GCA accession.
    Catches all errors so a bad accession doesn't abort the whole run.
    """
    row = {
        "gca": accession,
        "species_name": "",
        "taxon_id": "",
        "is_reference": "",
        "current_reference": "",
    }

    try:
        info = get_assembly_info(accession)
        row["taxon_id"] = info["taxon_id"]
        row["species_name"] = info["species_name"]

        ref_report = get_reference_for_taxon(str(info["taxon_id"]))

        if ref_report is None:
            row["is_reference"] = "No reference available"
            row["current_reference"] = "No reference available"
        else:
            ref_gca = ref_report.get("accession", "")
            ref_paired = ref_report.get("paired_accession", "")
            # current_reference: prefer the GCA accession
            row["current_reference"] = ref_gca or ref_paired or "Unknown"

            if accession in (ref_gca, ref_paired):
                row["is_reference"] = True
            else:
                row["is_reference"] = False

    except Exception as exc:  # noqa: BLE001
        row["is_reference"] = "Error"
        row["current_reference"] = f"Error: {exc}"
        print(f"  [error] {accession}: {exc}", file=sys.stderr)

    return row


def load_accessions_from_file(path: str) -> list[str]:
    """
    Accept a plain text file (one accession per line) or a CSV whose
    first column contains accessions (header row is skipped if it does
    not look like a GCA accession).
    """
    p = Path(path)
    if not p.exists():
        sys.exit(f"File not found: {path}")

    lines = p.read_text().splitlines()
    accessions = []
    for line in lines:
        # Strip BOM, whitespace, quotes
        val = line.strip().strip('"').strip("'")
        # Skip blank lines and header-like rows
        if not val or val.lower() in ("gca", "accession", "assembly"):
            continue
        # For CSV: take only the first column
        val = val.split(",")[0].strip()
        if val:
            accessions.append(val)
    return accessions


def write_csv(rows: list[dict], output_path: str) -> None:
    fieldnames = ["gca", "species_name", "taxon_id", "is_reference", "current_reference"]
    with open(output_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)



def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Check whether GCA accessions are reference assemblies (NCBI).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "accessions",
        nargs="*",
        metavar="GCA",
        help="One or more GCA accession(s) to check.",
    )
    parser.add_argument(
        "--file", "-f",
        metavar="PATH",
        help="Text or CSV file with one GCA accession per line / first column.",
    )
    parser.add_argument(
        "--output", "-o",
        default="reference_check.csv",
        metavar="PATH",
        help="Output CSV file path (default: reference_check.csv).",
    )
    return parser.parse_args()


def main() -> None:
    """Module entry point."""
    args = parse_args()

    accessions: list[str] = list(args.accessions)

    if args.file:
        accessions += load_accessions_from_file(args.file)

    if not accessions:
        sys.exit(
            "No accessions provided. "
            "Pass them on the command line or use --file."
        )

    # Deduplicate while preserving order
    seen: set[str] = set()
    unique: list[str] = []
    for a in accessions:
        if a not in seen:
            seen.add(a)
            unique.append(a)

    print(f"Checking {len(unique)} accession(s) …", file=sys.stderr)
    results = []
    for i, acc in enumerate(unique, 1):
        print(f"  [{i}/{len(unique)}] {acc}", file=sys.stderr)
        results.append(check_assembly(acc))

    write_csv(results, args.output)
    print(f"\nDone. Results written to: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()