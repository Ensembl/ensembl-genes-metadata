#!/usr/bin/env python3
"""Print ENA run metadata and linked EBI BioSamples records (no DB writes)."""
import argparse
import csv
import io
import json
import logging
import re
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen

PORTAL = 'https://www.ebi.ac.uk/ena/portal/api'
BROWSER = 'https://www.ebi.ac.uk/ena/browser/api/xml'
BIOSAMPLES = 'https://www.ebi.ac.uk/biosamples/samples'
LOG = logging.getLogger('ena_metadata')


def fetch(url, attempts=4):
    for attempt in range(attempts):
        try:
            req = Request(url, headers={'User-Agent': 'Ensembl-Genes-Metadata/1.0', 'Accept': '*/*'})
            with urlopen(req, timeout=90) as resp:
                return resp.read().decode('utf-8-sig')
        except (HTTPError, URLError, TimeoutError) as exc:
            if isinstance(exc, HTTPError) and exc.code not in (429, 500, 502, 503, 504):
                raise
            if attempt == attempts - 1:
                raise
            delay = min(2 ** attempt, 16)
            LOG.warning('Request failed (%s); retrying in %ss: %s', exc, delay, url)
            time.sleep(delay)


def portal_rows(params):
    url = PORTAL + '/search?' + urlencode(params)
    LOG.debug('ENA request: %s', url)
    body = fetch(url)
    rows = list(csv.DictReader(io.StringIO(body), delimiter='\t'))
    if rows and 'run_accession' not in rows[0]:
        raise ValueError('ENA response does not contain run_accession; response starts: ' + body[:250])
    return rows


def run_accessions_for_taxon(taxon):
    query = (f'tax_eq({taxon}) AND instrument_platform=ILLUMINA '
             'AND library_layout=PAIRED AND library_source=TRANSCRIPTOMIC')
    rows = portal_rows({'result': 'read_run', 'query': query, 'domain': 'read',
                        'fields': 'run_accession', 'limit': '0'})
    runs = sorted({r['run_accession'].strip() for r in rows if r.get('run_accession', '').strip()})
    LOG.info('Taxon %s: %d matching runs', taxon, len(runs))
    return runs


def available_fields():
    url = PORTAL + '/returnFields?' + urlencode({'result': 'read_run'})
    body = fetch(url)
    rows = list(csv.DictReader(io.StringIO(body), delimiter='\t'))
    fields = [r.get('columnId', '').strip() for r in rows if r.get('columnId')]
    if not fields:
        raise ValueError('Could not parse ENA returnFields response: ' + body[:250])
    return list(dict.fromkeys(['run_accession'] + fields))


def ena_run_metadata(run, fields):
    # ENA may reject very long requests: retry with core fields if necessary.
    params = {'result': 'read_run', 'query': f'run_accession="{run}"',
              'fields': ','.join(fields), 'format': 'tsv', 'limit': '0'}
    try:
        rows = portal_rows(params)
    except (HTTPError, ValueError) as exc:
        LOG.warning('Full field query failed for %s (%s); requesting core fields', run, exc)
        params['fields'] = ','.join(['run_accession', 'sample_accession',
            'secondary_sample_accession', 'study_accession', 'experiment_accession', 'experiment_alias','experiment_title',
            'tax_id', 'instrument_platform', 'library_strategy', 'library_layout','description',
            'library_name','tissue_type','tissue_lib','strain','cultivar','ecotype','isolate','cell_line','cell_type','fastq_ftp', 'fastq_md5', 'fastq_bytes'])
        rows = portal_rows(params)
    return next((r for r in rows if r.get('run_accession') == run), None)


def sample_ids_from_browser(run):
    """Find BioSamples (SAM*) identifiers from linked ENA XML, if available."""
    ids = []
    for accession in (run,):
        url = BROWSER + '/' + quote(accession, safe='')
        try:
            root = ET.fromstring(fetch(url))
        except (HTTPError, URLError, ET.ParseError) as exc:
            LOG.warning('ENA Browser XML unavailable for %s: %s', accession, exc)
            continue
        for elem in root.iter():
            for value in list(elem.attrib.values()) + ([elem.text] if elem.text else []):
                ids.extend(re.findall(r'\bSAM[END][A-Z0-9]+\b', str(value), flags=re.I))
            # ENA XML may include a sample accession (ERS/SRS/DRS) rather than SAM*.
            if elem.tag.upper().endswith('SAMPLE_REF'):
                ref = elem.get('accession') or elem.get('refname')
                if ref:
                    try:
                        sample_root = ET.fromstring(fetch(BROWSER + '/' + quote(ref, safe='')))
                        for e in sample_root.iter():
                            for value in list(e.attrib.values()) + ([e.text] if e.text else []):
                                ids.extend(re.findall(r'\bSAM[END][A-Z0-9]+\b', str(value), flags=re.I))
                    except (HTTPError, URLError, ET.ParseError) as exc:
                        LOG.warning('Cannot resolve ENA sample %s: %s', ref, exc)
    return list(dict.fromkeys(x.upper() for x in ids))


def biosample_for_run(metadata, run, cache):
    candidates = []
    for key in ('secondary_sample_accession', 'sample_accession'):
        for value in (metadata.get(key) or '').split(';'):
            value = value.strip()
            if re.fullmatch(r'SAM[END][A-Z0-9]+', value, flags=re.I):
                candidates.append(value.upper())
    if not candidates:
        candidates = sample_ids_from_browser(run)
    for accession in dict.fromkeys(candidates):
        if accession not in cache:
            try:
                cache[accession] = json.loads(fetch(BIOSAMPLES + '/' + quote(accession, safe='')))
            except (HTTPError, URLError, ValueError) as exc:
                LOG.warning('BioSamples lookup failed for %s: %s', accession, exc)
                cache[accession] = None
        if cache[accession] is not None:
            return accession, cache[accession]
    return (candidates[0] if candidates else None), None


def read_runs_file(path):
    p = Path(path).expanduser().resolve()
    if not p.is_file():
        raise FileNotFoundError(p)
    return [x.strip() for x in p.read_text(encoding='utf-8-sig').splitlines()
            if x.strip() and x.strip().lower() != 'run_accession']


def input_runs(args):
    runs = list(args.run or [])
    for taxon in args.taxon_id or []:
        runs.extend(run_accessions_for_taxon(taxon))
    if args.runs_file:
        runs.extend(read_runs_file(args.runs_file))
    if args.input_csv:
        csv_path = Path(args.input_csv).expanduser().resolve()
        with csv_path.open(newline='', encoding='utf-8-sig') as fh:
            reader = csv.DictReader(fh)
            if not reader.fieldnames or 'taxon_id' not in reader.fieldnames:
                raise ValueError('CSV must have taxon_id column')
            for row in reader:
                run_file = (row.get('runs_file') or row.get('run_file') or '').strip()
                if run_file:
                    p = Path(run_file).expanduser()
                    if not p.is_absolute():
                        p = csv_path.parent / p
                    runs.extend(read_runs_file(p))
                else:
                    taxon = (row.get('taxon_id') or '').strip()
                    if not taxon:
                        raise ValueError('CSV row missing taxon_id')
                    runs.extend(run_accessions_for_taxon(taxon))
    return list(dict.fromkeys(runs))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--taxon-id', action='append', help='Taxon ID (repeatable)')
    group.add_argument('--input-csv', help='small_taxa.csv or batch_001.csv')
    group.add_argument('--runs-file', help='Run accessions, one per line')
    group.add_argument('--run', action='append', help='Single run accession (repeatable)')
    parser.add_argument('--max-runs', type=int, default=0, help='Process at most N runs (0 = all)')
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format='%(levelname)s: %(message)s', stream=sys.stderr)
    try:
        runs = input_runs(args)
        if args.max_runs > 0:
            runs = runs[:args.max_runs]
        LOG.info('Processing %d distinct runs', len(runs))
        fields = available_fields()
        LOG.info('ENA supports %d indexed read_run fields', len(fields))
        cache = {}
        failures = 0
        for idx, run in enumerate(runs, 1):
            LOG.info('[%d/%d] %s', idx, len(runs), run)
            try:
                metadata = ena_run_metadata(run, fields)
                if metadata is None:
                    raise ValueError('Run not found in ENA Portal')
                biosample_accession, biosample = biosample_for_run(metadata, run, cache)
                record = {'run_accession': run, 'ena_indexed_metadata': metadata,
                          'biosample_accession': biosample_accession,
                          'biosample_metadata': biosample}
                print(json.dumps(record, ensure_ascii=False, indent=2), flush=True)
                if biosample is None:
                    LOG.warning('No BioSamples record resolved for %s', run)
            except (HTTPError, URLError, ValueError, TimeoutError) as exc:
                failures += 1
                LOG.error('Run %s failed: %s', run, exc)
        LOG.info('Finished: %d runs, %d failures, %d cached BioSamples', len(runs), failures, len(cache))
        return 1 if failures else 0
    except (OSError, ValueError, HTTPError, URLError) as exc:
        LOG.error('%s', exc)
        return 1


if __name__ == '__main__':
    sys.exit(main())
