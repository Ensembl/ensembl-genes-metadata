# ncbi_params.py

NCBI_PARAMS = {
    'bioprojects': None,
    'filters.reference_only': 'false',
    'filters.assembly_source': 'genbank',
    'filters.has_annotation': 'false',
    'filters.exclude_paired_reports': 'false',
    'filters.exclude_atypical': 'true',
    'filters.assembly_version': 'current',
    'filters.assembly_level': ['scaffold', 'chromosome', 'complete_genome', 'contig'],
    'filters.first_release_date': '04/04/2024',
    'filters.last_release_date': None,
    'filters.search_text': None,
    'filters.is_metagenome_derived': 'metagenome_derived_exclude',
    'tax_exact_match': None,
    'table_fields': ['assminfo-accession','assminfo-name'],
    'returned_content': 'ASSM_ACC',
    'page_size': '100', 
    'page_token': None, 
    'sort.field': None,
    'sort.direction': None,
    'include_tabular_header': 'INCLUDE_TABULAR_HEADER_ALWAYS'    
}