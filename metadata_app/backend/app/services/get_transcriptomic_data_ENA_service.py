import logging
import asyncio
import aiohttp
import pandas as pd

ENA_BASE_URL = "https://www.ebi.ac.uk/ena/portal/api/search?display=report&domain=read&result=read_run&fields=sample_accession,run_accession,fastq_ftp,read_count,instrument_platform"
MAX_RETRIES = 5  # Maximum number of retries for failed requests
CONCURRENT_LIMIT = 5  # Maximum concurrent API requests

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
                    elif response.status == 429:  # Too Many Requests (rate limited)
                        wait_time = 2 ** attempt  # Exponential backoff
                        print(f"Rate limited. Retrying in {wait_time} seconds...")
                        await asyncio.sleep(wait_time)
                    else:
                        response.raise_for_status()

        except (aiohttp.ClientResponseError, aiohttp.ClientConnectorError, aiohttp.ClientOSError) as e:
            wait_time = 2 ** attempt
            print(f"Request failed: {e}. Retrying in {wait_time} seconds...")
            await asyncio.sleep(wait_time)

    print(f"Failed to fetch data after {MAX_RETRIES} retries: {search_url}")
    return 0  # Return 0 if the request ultimately fails

async def check_data_from_ena(taxon_id, tree, semaphore):
    """Query ENA API for sequencing run counts using controlled concurrency."""
    query_base = f"tax_tree({taxon_id})" if tree else f"tax_eq({taxon_id})"

    queries = {
        "Short-read paired-end illumina": f"{query_base} AND instrument_platform=ILLUMINA AND library_layout=PAIRED AND library_source=TRANSCRIPTOMIC",
        "Short-read single-end illumina": f"{query_base} AND instrument_platform=ILLUMINA AND library_layout=SINGLE AND library_source=TRANSCRIPTOMIC",
        "Long-read PacBio": f"{query_base} AND instrument_platform=PACBIO_SMRT AND library_source=TRANSCRIPTOMIC",
        "Long-read ONP": f"{query_base} AND instrument_platform=OXFORD_NANOPORE AND library_source=TRANSCRIPTOMIC"
    }

    async with aiohttp.ClientSession() as session:
        tasks = {key: ena_rest_api(session, query, semaphore) for key, query in queries.items()}
        results = await asyncio.gather(*tasks.values())

    return {"Taxon ID": taxon_id, **dict(zip(queries.keys(), results))}

def add_data_from_ena(df):
    # Check transcriptomic data for each taxon_id in the dataset
    logging.info(f"Transcriptomic data check from ENA requested")
    # Get unique taxon IDs and filter out NaN values
    taxon_ids = [tid for tid in df["lowest_taxon_id"].unique() if pd.notna(tid)]
    species_taxon_ids = [tid for tid in df["species_taxon_id"].unique() if pd.notna(tid)]
    genus_taxon_ids = [gtid for gtid in df["genus_taxon_id"].unique() if pd.notna(gtid)]

    # Count NaN values and log warning if any are found
    nan_lowest_count = df["lowest_taxon_id"].isna().sum()
    nan_species_count = df["species_taxon_id"].isna().sum()
    nan_genus_count = df["genus_taxon_id"].isna().sum()

    if nan_lowest_count > 0:
        logging.warning(f"Found {nan_lowest_count} NA values in lowest_taxon_id column")

    if nan_lowest_count > 0:
        logging.warning(f"Found {nan_species_count} NA values in species_taxon_id column")

    if nan_genus_count > 0:
        logging.warning(f"Found {nan_genus_count} NA values in genus_taxon_id column")

    # Combine unique taxon IDs into a set and convert to integers
    all_taxon_ids = set()
    for tid in list(set(taxon_ids).union(set(genus_taxon_ids)).union(set(species_taxon_ids))):
        try:
            all_taxon_ids.add(int(tid))
        except (ValueError, TypeError) as e:
            logging.warning(f"Could not convert taxon ID {tid} to integer: {e}")

    if not all_taxon_ids:
        logging.warning("No valid taxon IDs found for transcriptomic data check")
    else:
        logging.info(f"Found {len(all_taxon_ids)} valid taxon IDs for transcriptomic data check")

    semaphore = asyncio.Semaphore(5)

    # Run transcriptomic data retrieval asynchronously
    async def fetch_transcriptomic_data():
        return await asyncio.gather(
            *[check_data_from_ena(taxon_id, tree=True, semaphore=semaphore) for taxon_id in all_taxon_ids]
        )

    transcriptomic_results = asyncio.run(fetch_transcriptomic_data())

    # Convert results into a DataFrame
    transcriptomic_df = pd.DataFrame(transcriptomic_results)

    logging.info(f"ENA check for transcriptomic data finished")

    return transcriptomic_df

