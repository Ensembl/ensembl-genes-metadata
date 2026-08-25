import asyncio
import json
import logging
import time
from pathlib import Path

import aiohttp
import pandas as pd

ENA_BASE_URL = "https://www.ebi.ac.uk/ena/portal/api/search?display=report&domain=read&result=read_run&fields=sample_accession,run_accession,fastq_ftp,read_count,instrument_platform"
MAX_RETRIES = 5
CONCURRENT_LIMIT = 5

BACKEND_ROOT = Path(__file__).resolve().parents[2]
CACHE_FILE = BACKEND_ROOT / "cache" / "ena_cache.json"
CACHE_TTL = 90 * 24 * 60 * 60


def load_cache():
    if CACHE_FILE.exists():
        try:
            with CACHE_FILE.open("r", encoding="utf-8") as handle:
                return json.load(handle)
        except json.JSONDecodeError:
            return {}
    return {}


def save_cache(cache):
    CACHE_FILE.parent.mkdir(parents=True, exist_ok=True)
    with CACHE_FILE.open("w", encoding="utf-8") as handle:
        json.dump(cache, handle)


async def ena_rest_api(session, query, semaphore):
    """Query the ENA API with rate limiting and retries."""
    search_url = f"{ENA_BASE_URL}&query={query}"

    for attempt in range(MAX_RETRIES):
        try:
            async with semaphore:
                async with session.get(search_url) as response:
                    if response.status == 200:
                        text = await response.text()
                        results = text.strip().split("\n")[1:]
                        return len(results)
                    if response.status == 429:
                        wait_time = 2 ** attempt
                        logging.warning(
                            "Rate limited by ENA. Retrying in %s seconds...",
                            wait_time,
                        )
                        await asyncio.sleep(wait_time)
                    else:
                        response.raise_for_status()

        except (
            aiohttp.ClientResponseError,
            aiohttp.ClientConnectorError,
            aiohttp.ClientOSError,
        ) as exc:
            wait_time = 2 ** attempt
            logging.warning(
                "ENA request failed: %s. Retrying in %s seconds...",
                exc,
                wait_time,
            )
            await asyncio.sleep(wait_time)

    logging.error("Failed to fetch ENA data after %s retries: %s", MAX_RETRIES, search_url)
    return 0


async def check_data_from_ena(taxon_id, tree, semaphore, cache, now):
    """Query ENA API for sequencing run counts using a shared cache."""
    cache_key = f"{taxon_id}:{tree}"

    if cache_key in cache and now - cache[cache_key]["timestamp"] < CACHE_TTL:
        return cache[cache_key]["data"]

    query_base = f"tax_tree({taxon_id})" if tree else f"tax_eq({taxon_id})"
    queries = {
        "short_read_paired_end_illumina": f"{query_base} AND instrument_platform=ILLUMINA AND library_layout=PAIRED AND library_source=TRANSCRIPTOMIC",
        "long_read_pacbio": f"{query_base} AND instrument_platform=PACBIO_SMRT AND library_source=TRANSCRIPTOMIC",
        "long_read_onp": f"{query_base} AND instrument_platform=OXFORD_NANOPORE AND library_source=TRANSCRIPTOMIC",
    }

    async with aiohttp.ClientSession() as session:
        tasks = {
            key: ena_rest_api(session, query, semaphore)
            for key, query in queries.items()
        }
        results = await asyncio.gather(*tasks.values())

    data = {"taxon_id": taxon_id, **dict(zip(queries.keys(), results))}
    cache[cache_key] = {"data": data, "timestamp": now}

    return data


def add_data_from_ena(df):
    """Check transcriptomic data for each taxon_id in the dataset."""
    logging.info("Transcriptomic data check from ENA requested")

    taxon_ids = {
        int(taxon_id)
        for taxon_id in pd.concat(
            [df["lowest_taxon_id"], df["species_taxon_id"], df["genus_taxon_id"]]
        )
        .dropna()
        .unique()
    }

    logging.info(
        "Found %s valid taxon IDs for transcriptomic data check",
        len(taxon_ids),
    )

    cache = load_cache()
    now = time.time()
    semaphore = asyncio.Semaphore(CONCURRENT_LIMIT)

    async def fetch_transcriptomic_data():
        tasks = [
            check_data_from_ena(
                taxon_id,
                tree=True,
                semaphore=semaphore,
                cache=cache,
                now=now,
            )
            for taxon_id in taxon_ids
        ]
        return await asyncio.gather(*tasks)

    transcriptomic_results = asyncio.run(fetch_transcriptomic_data())
    save_cache(cache)

    transcriptomic_df = pd.DataFrame(transcriptomic_results)
    logging.info("ENA check for transcriptomic data finished")
    return transcriptomic_df
