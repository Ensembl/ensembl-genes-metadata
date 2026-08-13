#!/usr/bin/env python3

import argparse
import json
import logging
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

API_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/dataset_report"
BACKEND_ROOT = Path(__file__).resolve().parents[2]
CACHE_FILE = BACKEND_ROOT / "cache" / "ncbi_reference_cache.json"
CACHE_TTL = 90 * 24 * 60 * 60  # 3 months in seconds


def load_cache():
    if CACHE_FILE.exists():
        try:
            with CACHE_FILE.open("r", encoding="utf-8") as handle:
                return json.load(handle)
        except json.JSONDecodeError:
            logging.warning("NCBI reference cache is invalid JSON. Rebuilding it.")
            return {}
    return {}


def save_cache(cache):
    CACHE_FILE.parent.mkdir(parents=True, exist_ok=True)
    with CACHE_FILE.open("w", encoding="utf-8") as handle:
        json.dump(cache, handle)


def get_ncbi_assembly(accession: str) -> dict:
    url = API_URL.format(accession=accession)

    request = urllib.request.Request(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": "ensembl-genes-metadata/1.0",
        },
    )

    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            return json.load(response)

    except urllib.error.HTTPError as exc:
        raise RuntimeError(f"NCBI API returned HTTP {exc.code} for {accession}") from exc

    except urllib.error.URLError as exc:
        raise RuntimeError(f"Could not connect to NCBI: {exc.reason}") from exc


def _build_reference_record(accession: str) -> dict:
    try:
        data = get_ncbi_assembly(accession)
        reports = data.get("reports", [])

        if not reports:
            raise RuntimeError(f"No assembly found for {accession}")

        report = reports[0]
        assembly_info = report.get("assembly_info", {})
        refseq_category = assembly_info.get("refseq_category")

        return {
            "gca": accession,
            "is_reference_genome": refseq_category == "reference genome",
            "refseq_category": refseq_category,
        }
    except RuntimeError as exc:
        logging.warning("NCBI reference lookup failed for %s: %s", accession, exc)
        return {
            "gca": accession,
            "is_reference_genome": None,
            "refseq_category": None,
        }


def get_reference_data_for_accessions(accessions: list[str]) -> pd.DataFrame:
    unique_accessions = sorted({accession for accession in accessions if accession})
    if not unique_accessions:
        return pd.DataFrame(columns=["gca", "is_reference_genome", "refseq_category"])

    cache = load_cache()
    now = time.time()
    records = []

    for accession in unique_accessions:
        cached_record = cache.get(accession)
        if cached_record and now - cached_record["timestamp"] < CACHE_TTL:
            records.append(cached_record["data"])
            continue

        record = _build_reference_record(accession)
        cache[accession] = {
            "data": record,
            "timestamp": now,
        }
        records.append(record)

    save_cache(cache)
    return pd.DataFrame(records)


def check_reference(accession: str):
    data = get_ncbi_assembly(accession)

    reports = data.get("reports", [])

    if not reports:
        raise RuntimeError(f"No assembly found for {accession}")

    report = reports[0]

    assembly_info = report.get("assembly_info", {})
    refseq_category = assembly_info.get("refseq_category")

    return refseq_category == "reference genome", refseq_category


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("accession")
    args = parser.parse_args()

    try:
        is_reference, category = check_reference(args.accession)

        print(f"Accession:       {args.accession}")
        print(f"RefSeq category: {category}")
        print(f"Is reference:    {is_reference}")

    except RuntimeError as exc:
        print(f"ERROR: {exc}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
