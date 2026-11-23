#!/usr/bin/env python3
"""Export ENA read_run and EBI BioSamples metadata into normalized CSVs; no DB writes."""
import argparse, csv, io, json, logging, re, sys, time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode, quote
from urllib.request import Request, urlopen

PORTAL='https://www.ebi.ac.uk/ena/portal/api'
BIOSAMPLES='https://www.ebi.ac.uk/biosamples/samples'
LOG=logging.getLogger('ena_export')
# Only fields needed by the four destination tables. Unsupported ENA fields
# are removed using /returnFields before querying.
ENA_FIELDS = [
    "run_accession", "sample_accession", "secondary_sample_accession", "study_accession",
    "fastq_ftp", "fastq_md5", "first_created", "instrument_platform",
    "tax_id", "library_strategy", "library_source", "library_selection",
    "library_layout", "library_name", "tissue_type", "sample_title",
    "sample_description", "experiment_title", "experiment_alias", "description",
    "scientific_name", "center_name", "study_title",
]

BIOSAMPLE_KEYS = [
    "accession", "taxId", "tissue", "age", "dev_stage", "sex",
    "organ", "sub_tissue", "title", "organism",
    "strain", "isolate", "cell_type", "cell_line", "cultivar",
    "ecotype", "disease", "infection_status", "treatment",
]

BIOSAMPLE_COLS=['sample_accession','tissue','age','dev_stage','sex','organ','sub_tissue','sample_title','experiment_title','description','sample_description','experiment_alias','scientific_name','strain','isolate','cell_type','cell_line','cultivar','ecotype']
CONDITION_COLS=['sample_accession','disease','infection_status','treatment']
RUN_COLS=['run_id','run_accession','sample_accession','study_accession','fastq_1','fastq_2','md5_1','md5_2','first_created','instrument_platform','taxon_id','library_strategy','library_source','library_selection','library_layout','library_name']
STUDY_COLS=['study_id','study_accession','center_name','study_title']

def fetch(url, attempts=5):
    for i in range(attempts):
        try:
            req=Request(url,headers={'User-Agent':'RNAseq-Metadata-Exporter/1.0','Accept':'application/json, text/tab-separated-values, */*'})
            with urlopen(req,timeout=120) as r: return r.read().decode('utf-8-sig')
        except (HTTPError,URLError,TimeoutError) as e:
            if isinstance(e,HTTPError) and e.code not in (429,500,502,503,504): raise
            if i==attempts-1: raise
            wait=min(2**i,30)
            LOG.warning('Retry %s after %ss: %s',url,wait,e)
            time.sleep(wait)

def portal(params):
    body=fetch(PORTAL+'/search?'+urlencode(params))
    if body.lstrip().startswith(('{','[')):
        raise ValueError('Unexpected ENA JSON/error response: '+body[:200])
    return csv.DictReader(io.StringIO(body),delimiter='\t')

def field_names():
    rows=csv.DictReader(io.StringIO(fetch(PORTAL+'/returnFields?result=read_run')),delimiter='\t')
    return {r.get('columnId','') for r in rows}

def clean(v):
    if v is None: return ''
    if isinstance(v,list): return '; '.join(filter(None,(clean(x) for x in v)))
    if isinstance(v,dict): return clean(v.get('text') or v.get('value') or '')
    return str(v).strip()

def characteristics(bs):
    """Retain BioSamples characteristic keys exactly as returned."""
    result = {}
    for key, items in (bs.get('characteristics') or {}).items():
        vals = items if isinstance(items, list) else [items]
        value = '; '.join(filter(None, (clean(x) for x in vals)))
        if value:
            result[key] = value
    return result

MISSING_MARKERS = {'', 'missing', 'unknown', 'not determined', 'not provided',
                   'not available', 'not applicable', 'unspecified', 'na', 'n/a',
                   'null', 'none', 'no data', 'not collected', 'not recorded'}

def usable(value):
    return bool(clean(value)) and clean(value).casefold() not in MISSING_MARKERS

def choose(chars, ena, field):
    """Prefer usable exact-name BioSamples characteristic, then exact ENA field."""
    b = clean(chars.get(field))
    if field == 'tissue':
        return b  # BioSamples tissue only, never ENA tissue_type.
    e = clean(ena.get(field))
    if usable(b):
        return b
    if usable(e):
        return e
    return b or e



def biosample_value(bs, chars, ena, field):
    """Only explicitly specified source-to-output mappings; no fuzzy aliases."""

    # Tissue: BioSamples only
    if field == 'tissue':
        return clean(chars.get('tissue'))

    # Description: ENA only
    if field == 'description':
        return clean(ena.get('description'))

    # Experiment metadata: ENA only
    if field in ('experiment_title', 'experiment_alias'):
        return clean(ena.get(field))

    # Sample title: BioSamples first, ENA fallback
    if field == 'sample_title':
        b = clean(chars.get('title'))
        e = clean(ena.get('sample_title'))

    # Sample description: BioSamples first, ENA fallback
    elif field == 'sample_description':
        b = clean(chars.get('description'))
        e = clean(ena.get('sample_description'))

    # Scientific name: BioSamples first, ENA fallback
    elif field == 'scientific_name':
        b = clean(chars.get('organism'))
        e = clean(ena.get('scientific_name'))

    # All other fields: exact-name matching
    else:
        return choose(chars, ena, field)

    if usable(b):
        return b
    if usable(e):
        return e

    return b or e


def get_biosample(ena,cache):
    candidates=[]
    for key in ('secondary_sample_accession','sample_accession'):
        for value in clean(ena.get(key)).split(';'):
            value=value.strip().upper()
            if re.fullmatch(r'SAM[END][A-Z0-9]+',value): candidates.append(value)
    for accession in dict.fromkeys(candidates):
        if accession not in cache:
            try: cache[accession]=json.loads(fetch(BIOSAMPLES+'/'+quote(accession,safe='')))
            except (HTTPError,URLError,ValueError) as e:
                LOG.warning('BioSamples lookup failed %s: %s',accession,e)
                cache[accession]={}
        if cache[accession]: return cache[accession]
    return {}

def split_fastq(ena):
    urls=[x.strip() for x in clean(ena.get('fastq_ftp')).split(';') if x.strip()]
    md5s=[x.strip() for x in clean(ena.get('fastq_md5')).split(';')]
    if len(urls)>2 or len(md5s)>2: raise ValueError('More than two FASTQ/MD5 values for run %s; refusing to discard files' % ena.get('run_accession'))
    if urls and len(md5s)!=len(urls): raise ValueError('FASTQ/MD5 count mismatch: %s' % ena.get('run_accession'))
    return dict(fastq_1=urls[0] if urls else '',fastq_2=urls[1] if len(urls)>1 else '',md5_1=md5s[0] if md5s else '',md5_2=md5s[1] if len(md5s)>1 else '')

def runs_from_file(path):
    return list(dict.fromkeys(x.strip() for x in Path(path).read_text(encoding='utf-8-sig').splitlines() if x.strip() and not x.startswith('#') and x.strip().lower()!='run_accession'))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--taxon-id',action='append',help='NCBI taxon ID; repeatable')
    g.add_argument('--runs-file',help='Run accessions, one per line')
    g.add_argument('--run',action='append',help='Run accession; repeatable')
    p.add_argument('--outdir',default='metadata_export')
    p.add_argument('--batch-size',type=int,default=1000,help='ENA query page size (default 1000)')
    p.add_argument('--max-runs',type=int,default=0)
    p.add_argument('--skip-existing',action='store_true',help='Append only unseen runs from existing run_metadata.csv')
    p.add_argument('--verbose',action='store_true')
    a=p.parse_args()
    logging.basicConfig(level=logging.DEBUG if a.verbose else logging.INFO,format='%(levelname)s: %(message)s',stream=sys.stderr)
    out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True)
    available=field_names()
    fields=[f for f in ENA_FIELDS if f in available]
    required={'run_accession','sample_accession','study_accession'}
    if not required.issubset(fields): raise ValueError('ENA missing required fields: '+str(required-set(fields)))
    LOG.info('Requesting %d supported ENA fields; skipped unsupported: %s',len(fields), ', '.join(f for f in ENA_FIELDS if f not in available) or 'none')
    seen_runs=set();seen_samples=set();seen_conditions=set();seen_studies=set()
    if a.skip_existing:
        for name,key,target in [('run_metadata.csv','run_accession',seen_runs),('biosample_metadata.csv','sample_accession',seen_samples),('sample_condition_metadata.csv','sample_accession',seen_conditions),('study_metadata.csv','study_accession',seen_studies)]:
            path=out/name
            if path.exists():
                with path.open(newline='',encoding='utf-8') as f:
                    target.update(r[key] for r in csv.DictReader(f) if r.get(key))
    outputs={}
    for name,columns in [('biosample',BIOSAMPLE_COLS),('sample_condition',CONDITION_COLS),('run',RUN_COLS),('study',STUDY_COLS)]:
        path=out/(name+'_metadata.csv')
        append=a.skip_existing and path.exists() and path.stat().st_size>0
        fh=path.open('a' if append else 'w',newline='',encoding='utf-8')
        writer=csv.DictWriter(fh,fieldnames=columns,extrasaction='ignore');
        if not append: writer.writeheader()
        outputs[name]=(fh,writer)
    cache={};count=0
    try:
        if a.taxon_id:
            def source():
                for taxon in a.taxon_id:
                    query=f'tax_eq({taxon}) AND instrument_platform=ILLUMINA AND library_layout=PAIRED AND library_source=TRANSCRIPTOMIC'
                    offset=0
                    while True:
                        rows=list(portal({'result':'read_run','query':query,'fields':','.join(fields),'limit':str(a.batch_size),'offset':str(offset),'format':'tsv'}))
                        if not rows: break
                        yield from rows
                        offset+=len(rows)
                        LOG.info('Taxon %s: fetched %d ENA records',taxon,offset)
                        if len(rows)<a.batch_size: break
            records=source()
        else:
            ids=a.run or runs_from_file(a.runs_file)
            def source_ids():
                for run in dict.fromkeys(ids):
                    rows=list(portal({'result':'read_run','query':f'run_accession="{run}"','fields':','.join(fields),'format':'tsv'}))
                    if not rows: LOG.warning('Run not found: %s',run)
                    yield from rows
            records=source_ids()
        for ena in records:
            run=clean(ena.get('run_accession'))
            if not run or run in seen_runs: continue
            sample=clean(ena.get('sample_accession'))
            study=clean(ena.get('study_accession'))
            bs=get_biosample(ena,cache)
            # sample_accession remains ENA's sample ID (ERS/SRS/DRS).
            # BioSamples SAM accession is not interchangeable with that FK.
            chars=characteristics(bs)
            if sample and sample not in seen_samples:
                row={'sample_accession':sample}
                for key in BIOSAMPLE_COLS[1:]: row[key]=biosample_value(bs,chars,ena,key)
                outputs['biosample'][1].writerow(row);seen_samples.add(sample)
            if sample and sample not in seen_conditions:
                row={'sample_accession':sample}
                for key in CONDITION_COLS[1:]: row[key]=choose(chars,ena,key)
                outputs['sample_condition'][1].writerow(row);seen_conditions.add(sample)
            if study and study not in seen_studies:
                outputs['study'][1].writerow({'study_accession':study,'center_name':clean(ena.get('center_name')),'study_title':clean(ena.get('study_title'))})
                seen_studies.add(study)
            row={key:clean(ena.get(key)) for key in RUN_COLS}
            # tax_id comes from ENA only; do not rename BioSamples taxId.
            # Map ENA tax_id to database taxonomy_id
            row['taxon_id'] = clean(ena.get('tax_id'))
            row.update(split_fastq(ena));outputs['run'][1].writerow(row)
            seen_runs.add(run);count+=1
            if count%100==0:
                for fh,_ in outputs.values():fh.flush()
                LOG.info('Exported %d new runs (%d BioSamples cached)',count,len(cache))
            if a.max_runs and count>=a.max_runs:break
    finally:
        for fh,_ in outputs.values():fh.close()
    LOG.info('Done: %d new runs; %d distinct samples; %d studies',count,len(seen_samples),len(seen_studies))

if __name__=='__main__':
    try:main()
    except (HTTPError,URLError,ValueError,OSError) as exc:
        LOG.error('%s',exc);sys.exit(1)
