# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

""" Script to check for the availability of transcriptomic data when provided with the taxon_id 
    This would also download meta data anf format the required fields into a tab delimited file """

import sys
import argparse
import re
import requests

""" Use taxon_id to check for the existence of short reads (paired ended) or long reads from the ENA """
def transcriptomic_status(taxon_id,output_path):
    """ Check status of transcriptomic data both in the registry for exisiting species and publicly for new species """
    ena_base_url = 'https://www.ebi.ac.uk/ena/portal/api/search?display=report'
    files_domain = 'domain=read&result=read_run'
    sample_domain = 'domain=sample&result=sample'
    sample_fields = 'accession,secondary_sample_accession,bio_material,cell_line,cell_type,collected_by,collection_date,country,cultivar,culture_collection,description,dev_stage,ecotype,environmental_sample,first_public,germline,identified_by,isolate,isolation_source,location,mating_type,serotype,serovar,sex,submitted_sex,specimen_voucher,strain,sub_species,sub_strain,tissue_lib,tissue_type,variety,tax_id,scientific_name,sample_alias,center_name,protocol_label,project_name,investigation_type,experimental_factor,sample_collection,sequencing_method'
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
                    sample = read_entry.split('\t')
                    #Query to get further info about the sample
                    print('Looking for meta data for sample '+sample[3])
                    sample_url = ena_base_url + 'query="accession=' + sample[3] +'"' + sample_domain + 'fields=' + sample_fields
                    print('Url is '+sample_url)
                    response = requests.get(sample_url)
                    encoding_style = response.encoding
                    if response:
                        print('We found content')
                        sample_meta = output_path + "/samples/" + sample[3]
                        open(sample_meta, "wb").write(response.content)
        if entry > 0:
            status = 'available'
    return status

""" Main method to handle arguments """
if __name__ == '__main__':
    """Main script entry-point."""
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_path', help='Path to store transcriptomic data report', required=True)
    parser.add_argument('--taxon_id', help='Taxonomy id to fetch information about', required=True)
    args = parser.parse_args()
    taxon_id = args.taxon_id
    output_path = args.output_path
    status = transcriptomic_status(taxon_id,output_path)
    print("Transcriptomic data  is "+str(status)+'\n')
