""" Script to check for the availability of transcriptomic data when provided with the taxon_id """

import sys
import argparse
import re
import requests

def transcriptomic_status(taxon_id,output_path):
    """ Check status of transcriptomic data both in the registry for exisiting species and publicly for new species """
    ena_base_url = 'https://www.ebi.ac.uk/ena/portal/api/search?display=report'
    files_domain = 'domain=read&result=read_run'
    files_fields = 'run_accession,study_accession,experiment_accession,sample_accession,secondary_sample_accession,instrument_platform,instrument_model,library_layout,library_strategy,nominal_length,read_count,base_count,fastq_ftp,fastq_aspera,submitted_ftp,fastq_md5,library_source,library_selection,center_name,study_alias,experiment_alias,experiment_title,study_title'
    download_method = 'ftp'
    total_read = 0
    entry = 0
    status = 'not available'
    ftp_path = ena_base_url + '&query="tax_tree(' + taxon_id + ') AND instrument_platform=ILLUMINA AND library_source=TRANSCRIPTOMIC"&' + files_domain + '&fields=' + files_fields
    response = requests.get(ftp_path)
    fastq_file = 'fastq_' + download_method
    header_pattern = '^(\w+.*)$'
    content_pattern = '^[a-z]'
    encoding_style = response.encoding
    unwanted_reads = []
    wanted_reads = []

    if response:
        read_file = output_path + "/ena_read_set.txt"
        open(read_file, "wb").write(response.content)
        with open(read_file) as f:
            for line in f.readlines():
                header = re.search(header_pattern, line)
                if header:
                    read_entry = header.group(1)
                content = re.search(content_pattern, read_entry)
                if not content:
                    total_read += 1
                    wanted_reads.append(read_entry)
                    entry += 1
        if entry > 0:
            status = 'available'
    return status


if __name__ == '__main__':
    """Main script entry-point."""
    parser = argparse.ArgumentParser()
    parser.add_argument('-out','--output_path', help='Path to store transcriptomic data report', required=True)
    parser.add_argument('-tid','--taxon_id', help='Taxonomy id to fetch information about', required=True)
    args = parser.parse_args()
    taxon_id = args.taxon_id
    output_path = args.output_path
    status = transcriptomic_status(taxon_id,output_path)
    print("Transcriptomic data  is "+str(status)+'\n')
