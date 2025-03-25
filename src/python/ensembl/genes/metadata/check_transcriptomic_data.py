import argparse
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

async def main():
    """Entrypoint"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--taxon_ids", type=str, nargs="+", help="List of taxon IDs (overrides CSV)")
    parser.add_argument("-f", "--file", type=str, help="CSV file containing a 'lowest_taxon_id' column")
    parser.add_argument("-o", "--output", type=str, help="Output CSV file to save results")
    parser.add_argument("--tree", action="store_true", help="Include descendant taxa in ENA query")
    args = parser.parse_args()

    # Read taxon IDs from CSV file if provided
    taxon_ids = []
    if args.file:
        df = pd.read_csv(args.file)
        if "lowest_taxon_id" not in df.columns:
            raise ValueError("CSV file must contain a column named 'lowest_taxon_id'")
        taxon_ids = df["lowest_taxon_id"].dropna().astype(str).tolist()

    # Override with CLI taxon IDs if provided
    if args.taxon_ids:
        taxon_ids = args.taxon_ids

    if not taxon_ids:
        raise ValueError("No taxon IDs provided. Use -f <file.csv> or -t <taxon_id>")

    # Limit concurrent API calls using Semaphore
    semaphore = asyncio.Semaphore(CONCURRENT_LIMIT)

    # Run all taxon ID queries with controlled concurrency
    results = await asyncio.gather(*[check_data_from_ena(taxon_id, args.tree, semaphore) for taxon_id in taxon_ids])

    # Convert results to a Pandas DataFrame
    df = pd.DataFrame(results)
    print(df)  # Print to console

    # Save to CSV if an output file is provided
    if args.output:
        df.to_csv(args.output, index=False)
        print(f"Results saved to {args.output}")

# Run script in an async event loop
if __name__ == "__main__":
    import sys
    if sys.argv[1:]:  # Only run CLI mode if arguments are provided
        asyncio.run(main())