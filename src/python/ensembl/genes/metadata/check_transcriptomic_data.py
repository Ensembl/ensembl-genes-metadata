import argparse
import asyncio
import aiohttp
import pandas as pd

ENA_BASE_URL = "https://www.ebi.ac.uk/ena/portal/api/search?display=report&domain=read&result=read_run&fields=sample_accession,run_accession,fastq_ftp,read_count,instrument_platform"
MAX_RETRIES = 5  # Maximum number of retries for failed requests

async def ena_rest_api(session, query):
    """Asynchronous function to query ENA API and count sequencing runs."""
    search_url = f"{ENA_BASE_URL}&query={query}"

    for attempt in range(MAX_RETRIES):
        try:
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

async def check_data_from_ena(taxon_id, tree):
    """Query ENA API for sequencing run counts using asynchronous requests."""
    query_base = f"tax_tree({taxon_id})" if tree else f"tax_eq({taxon_id})"

    queries = {
        "Short-read paired-end illumina": f"{query_base} AND instrument_platform=ILLUMINA AND library_layout=PAIRED AND library_source=TRANSCRIPTOMIC",
        "Short-read single-end illumina": f"{query_base} AND instrument_platform=ILLUMINA AND library_layout=SINGLE AND library_source=TRANSCRIPTOMIC",
        "Long-read PacBio": f"{query_base} AND instrument_platform=PACBIO_SMRT AND library_source=TRANSCRIPTOMIC",
        "Long-read ONP": f"{query_base} AND instrument_platform=OXFORD_NANOPORE AND library_source=TRANSCRIPTOMIC"
    }

    async with aiohttp.ClientSession() as session:
        tasks = {key: ena_rest_api(session, query) for key, query in queries.items()}
        results = await asyncio.gather(*tasks.values())

    return {"Taxon ID": taxon_id, **dict(zip(queries.keys(), results))}

async def main():
    """Entrypoint"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--taxon_ids", type=str, nargs="+", required=True, help="List of taxon IDs")
    parser.add_argument("--tree", action="store_true", help="Include descendant taxa in ENA query")
    args = parser.parse_args()

    # Run all taxon ID queries concurrently
    results = await asyncio.gather(*[check_data_from_ena(taxon_id, args.tree) for taxon_id in args.taxon_ids])

    # Convert results to a Pandas DataFrame and display
    df = pd.DataFrame(results)
    print(df)

# Run script in an async event loop
if __name__ == "__main__":
    import sys
    if sys.argv[1:]:  # Only run CLI mode if arguments are provided
        asyncio.run(main())