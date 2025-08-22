import logging
import asyncio
import aiohttp
import pandas as pd
import os
import json
import time

ENA_BASE_URL = "https://www.ebi.ac.uk/ena/portal/api/search?display=report&domain=read&result=read_run&fields=sample_accession,run_accession,fastq_ftp,read_count,instrument_platform"
MAX_RETRIES = 5  # Maximum number of retries for failed requests
CONCURRENT_LIMIT = 5  # Maximum concurrent API requests

CACHE_FILE = "metadata_app/backend/cache/ena_cache.json"
CACHE_TTL = 90 * 24 * 60 * 60  #3 months in seconds


def load_cache():
    if os.path.exists(CACHE_FILE):
        try:
            with open(CACHE_FILE, "r") as f:
                return json.load(f)
        except json.JSONDecodeError:
            return {}
    return {}


def save_cache(cache):
    with open(CACHE_FILE, "w") as f:
        json.dump(cache, f)


async def ena_rest_api(session, query, semaphore):
    """Asynchronous function to query ENA API with rate limiting and retries."""
    search_url = f"{ENA_BASE_URL}&query={query}"

    for attempt in range(MAX_RETRIES):
        try:
            async with semaphore:  # Limit concurrent tasks
                async with session.get(search_url) as response:
                    if response.status == 200:
                        text = await response.text()
                        results = text.strip().split("\n")[1:]  # Ignore header row
                        return len(results)
                    elif response.status == 429:  # Too Many Requests
                        wait_time = 2 ** attempt
                        logging.warning(f"Rate limited. Retrying in {wait_time} seconds...")
                        await asyncio.sleep(wait_time)
                    else:
                        response.raise_for_status()

        except (aiohttp.ClientResponseError, aiohttp.ClientConnectorError, aiohttp.ClientOSError) as e:
            wait_time = 2 ** attempt
            logging.warning(f"Request failed: {e}. Retrying in {wait_time} seconds...")
            await asyncio.sleep(wait_time)

    logging.error(f"Failed to fetch data after {MAX_RETRIES} retries: {search_url}")
    return 0


async def check_data_from_ena(taxon_id, tree, semaphore, cache, now):
    """Query ENA API for sequencing run counts using controlled concurrency + shared cache dict."""
    cache_key = f"{taxon_id}:{tree}"

    if cache_key in cache and now - cache[cache_key]["timestamp"] < CACHE_TTL:
        return cache[cache_key]["data"]

    query_base = f"tax_tree({taxon_id})" if tree else f"tax_eq({taxon_id})"
    queries = {
        "Short-read paired-end illumina": f"{query_base} AND instrument_platform=ILLUMINA AND library_layout=PAIRED AND library_source=TRANSCRIPTOMIC",
        "Long-read PacBio": f"{query_base} AND instrument_platform=PACBIO_SMRT AND library_source=TRANSCRIPTOMIC",
        "Long-read ONP": f"{query_base} AND instrument_platform=OXFORD_NANOPORE AND library_source=TRANSCRIPTOMIC"
    }

    async with aiohttp.ClientSession() as session:
        tasks = {key: ena_rest_api(session, query, semaphore) for key, query in queries.items()}
        results = await asyncio.gather(*tasks.values())

    data = {"Taxon ID": taxon_id, **dict(zip(queries.keys(), results))}
    cache[cache_key] = {"data": data, "timestamp": now}

    return data


def add_data_from_ena(df):
    """Check transcriptomic data for each taxon_id in the dataset (cached)."""
    logging.info("Transcriptomic data check from ENA requested")

    taxon_ids = [tid for tid in df["lowest_taxon_id"].unique() if pd.notna(tid)]
    species_taxon_ids = [tid for tid in df["species_taxon_id"].unique() if pd.notna(tid)]
    genus_taxon_ids = [gtid for gtid in df["genus_taxon_id"].unique() if pd.notna(gtid)]

    nan_lowest_count = df["lowest_taxon_id"].isna().sum()
    nan_species_count = df["species_taxon_id"].isna().sum()
    nan_genus_count = df["genus_taxon_id"].isna().sum()

    if nan_lowest_count > 0:
        logging.warning(f"Found {nan_lowest_count} NA values in lowest_taxon_id column")
    if nan_species_count > 0:
        logging.warning(f"Found {nan_species_count} NA values in species_taxon_id column")
    if nan_genus_count > 0:
        logging.warning(f"Found {nan_genus_count} NA values in genus_taxon_id column")

    all_taxon_ids = {int(tid) for tid in set(taxon_ids) | set(species_taxon_ids) | set(genus_taxon_ids) if
                     pd.notna(tid)}
    logging.info(f"Found {len(all_taxon_ids)} valid taxon IDs for transcriptomic data check")

    cache = load_cache()
    now = time.time()
    semaphore = asyncio.Semaphore(CONCURRENT_LIMIT)

    async def fetch_transcriptomic_data():
        return await asyncio.gather(
            *[check_data_from_ena(taxon_id, tree=True, semaphore=semaphore, cache=cache, now=now)
              for taxon_id in all_taxon_ids]
        )

    transcriptomic_results = asyncio.run(fetch_transcriptomic_data())

    # ✅ Save cache only once after all requests
    save_cache(cache)

    transcriptomic_df = pd.DataFrame(transcriptomic_results)
    logging.info("ENA check for transcriptomic data finished")
    return transcriptomic_df