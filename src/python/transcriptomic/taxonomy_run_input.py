#!/usr/bin/env python3
"""Retrieve ENA Illumina paired transcriptomic run accessions and batch large taxa."""
import argparse
import csv
import logging
import re
import sys
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

ENA_URL = 'https://www.ebi.ac.uk/ena/portal/api/search'
LOG = logging.getLogger('ena-batches')
RUN_PATTERN = re.compile(r'^[SED]RR\d+$')


def fetch_runs(taxon_id, retries=4, timeout=120):
    query = (f'tax_eq({taxon_id}) AND instrument_platform=ILLUMINA '
             'AND library_layout=PAIRED AND library_source=TRANSCRIPTOMIC')
    params = {'result': 'read_run', 'query': query, 'domain': 'read',
              'fields': 'run_accession', 'limit': '0'}
    url = f'{ENA_URL}?{urlencode(params)}'
    LOG.debug('ENA request for taxon %s: %s', taxon_id, url)
    for attempt in range(retries + 1):
        try:
            request = Request(url, headers={'User-Agent': 'Ensembl-Genes-Metadata/1.0',
                                            'Accept': 'text/plain'})
            with urlopen(request, timeout=timeout) as response:
                body = response.read().decode('utf-8-sig')
            lines = [line.strip() for line in body.splitlines() if line.strip()]
            if not lines or lines[0].lower() != 'run_accession':
                raise ValueError(f'Unexpected ENA response for taxon {taxon_id}: {lines[:3]!r}')
            runs = set(lines[1:])
            invalid = sorted(run for run in runs if not RUN_PATTERN.fullmatch(run))
            if invalid:
                raise ValueError(f'Unexpected run accession(s) for taxon {taxon_id}: {invalid[:5]}')
            return sorted(runs)
        except (HTTPError, URLError, TimeoutError) as exc:
            retryable = not isinstance(exc, HTTPError) or exc.code in (429, 500, 502, 503, 504)
            if not retryable or attempt == retries:
                raise RuntimeError(f'ENA request failed for taxon {taxon_id}: {exc}') from exc
            delay = min(2 ** attempt, 30)
            LOG.warning('ENA taxon %s request failed (%s); retrying in %ss', taxon_id, exc, delay)
            time.sleep(delay)
    raise AssertionError('Unreachable')


def write_csv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def load_taxa(path):
    with path.open(newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or not {'taxon_id', 'gca'}.issubset(reader.fieldnames):
            raise ValueError('Input CSV must contain taxon_id,gca columns')
        pairs = []
        seen = {}
        for line_number, row in enumerate(reader, start=2):
            taxon = (row.get('taxon_id') or '').strip()
            gca = (row.get('gca') or '').strip()
            if not taxon or not gca or not taxon.isdigit():
                raise ValueError(f'Invalid taxon_id/gca at CSV line {line_number}')
            if taxon in seen and seen[taxon] != gca:
                raise ValueError(f'Taxon {taxon} has multiple GCAs ({seen[taxon]}, {gca}); taxon-only folders would collide')
            if taxon not in seen:
                pairs.append((taxon, gca))
                seen[taxon] = gca
    return pairs


def run(input_csv, outdir, batch_size, retries, timeout):
    if batch_size < 1:
        raise ValueError('--batch-size must be greater than zero')
    output = outdir.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    taxa = load_taxa(input_csv)
    small = []
    counts = []
    for taxon, gca in taxa:
        runs = fetch_runs(taxon, retries=retries, timeout=timeout)
        count = len(runs)
        batch_count = 0 if count <= batch_size else (count + batch_size - 1) // batch_size
        counts.append({'taxon_id': taxon, 'gca': gca, 'run_count': count, 'batches': batch_count})
        if count <= batch_size:
            small.append({'taxon_id': taxon, 'gca': gca})
            print(f'{taxon} ({gca}): {count} runs -> small_taxa.csv')
            continue
        folder = output / 'batches' / taxon
        folder.mkdir(parents=True, exist_ok=True)
        for batch_number, start in enumerate(range(0, count, batch_size), start=1):
            subset = runs[start:start + batch_size]
            stem = f'batch_{batch_number:03d}'
            runs_path = folder / f'{stem}_runs.txt'
            runs_path.write_text('\n'.join(subset) + '\n', encoding='utf-8')
            write_csv(folder / f'{stem}.csv', ['taxon_id', 'gca', 'runs_file'],
                      [{'taxon_id': taxon, 'gca': gca, 'runs_file': str(runs_path.resolve())}])
        print(f'{taxon} ({gca}): {count} runs -> {batch_count} batches in {folder}')
    write_csv(output / 'small_taxa.csv', ['taxon_id', 'gca'], small)
    write_csv(output / 'run_counts.csv', ['taxon_id', 'gca', 'run_count', 'batches'], counts)
    print(f'\nDone. Output directory: {output}')
    print(f'Small taxa: {len(small)}; large taxa: {sum(row["batches"] > 0 for row in counts)}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True, help='Input CSV with taxon_id,gca')
    parser.add_argument('--outdir', type=Path, required=True, help='Output directory (absolute paths written to batch CSVs)')
    parser.add_argument('--batch-size', type=int, default=500)
    parser.add_argument('--retries', type=int, default=4)
    parser.add_argument('--timeout', type=int, default=120)
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format='%(levelname)s: %(message)s')
    try:
        run(args.input, args.outdir, args.batch_size, args.retries, args.timeout)
    except (OSError, ValueError, RuntimeError) as exc:
        LOG.error('%s', exc)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
