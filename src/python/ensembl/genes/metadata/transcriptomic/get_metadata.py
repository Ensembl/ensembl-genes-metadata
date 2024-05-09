#  See the NOTICE file distributed with this work for additional information
#  regarding copyright ownership.
#
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#  http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
"""Script to produce a json file per given run, containing all accessions, in the format:

{
    table : [
        metakey : metavalue,
        ...
       ]
},...

Returns:
    str : json element containing metadata for given run
"""

import argparse
import requests
import json

def _request_data(run_accession,fields) -> str:
    """This function, makes a HTTP request for the metadata of the given run_accession.

    Args:
        run_accession (str): unique run aceession
        fields (array): list of fields of interest to have

    Returns:
        str: tsv file containing said metadata
    """

    query = 'run_accession="' + run_accession + '"'

    data = {
    'result': 'read_run',
    'query': query,
    'fields': ",".join(fields),
    'format': 'tsv',
    }

    request = requests.post('https://www.ebi.ac.uk/ena/portal/api/search', data=data)
    response = request.text

    return response

def _json_parse(response,fields) -> str:
    """Parse response from http request into a readable json

    Args:
        data (str): tsv from the http request
        fields (array): list of fields of interest to have

    Returns:
        str: json formatted dictionary of metadata and values
    """
    
    table = response.split('\n')
    data = table[1].split('\t')

    d = {}
    for key in fields: 
        for value in data:
            d[key] = value
            data.remove(value)
            break

    # we load the data from before in an orderly fashion to fit out format
    table_run = {
        # run_id -> vianey's script get it and assign it
        'taxon_id' : d['tax_id'],
        'run_accession' : d['run_accession'],
        # species (when subspecies) and genus tax_id -> not in the list of available fields...
        'qc_status' : 'NOT_CHECKED',
        'sample_accession' : d['sample_accession'],
        'study_accession' : d['study_accession'],
        'read_type' : d['library_strategy'],
        'platform' : d['instrument_platform'],
        'paired' : 1 if d['library_layout'] == "PAIRED" else 0,
        'experiment' : '; '.join(value for value in [d['experiment_alias'],d['experiment_title']] if value is not None),
        'description' : d['description'],
        'library_name' : '; '.join(value for value in [d['library_source'],d['library_name']] if value is not None),
        'library_selection' : d['library_selection'],
        'tissue' : '; '.join(value for value in [d['tissue_type'],d['tissue_lib']] if value != ""),
        'cell_line' : d['cell_line'],
        'cell_type' : d['cell_type'],
        'strain' : '; '.join(value for value in [d['strain'],d['cultivar'],d['ecotype'],d['isolate']] if value != ""),
    }

    table_study = {
        'study_accession' : d['study_accession'],
        'center_name' : d['center_name'],
    }

    table_data_files = list()

    file_name = list()
    file_url = d['fastq_ftp'].split(';')
    file_md5 = d['fastq_md5'].split(';')

    for url in (file_url):
        extension_name = url.split('/')[-1]
        base_name = extension_name.split('.')[0]
        file_name.append(base_name)

    for i in range(len(file_name)):
        read = {
            'name'  : file_name[i],
            'url'   : file_url[i],
            'md5'   : file_md5[i]
        }
        table_data_files.append(read)

    json_run = json.dumps(table_run)
    json_study = json.dumps(table_study)
    json_data_files = json.dumps(table_data_files)

    with open("insert_into_run.json", "w") as file:
        file.write(json_run)
    with open("insert_into_study.json", "w") as file:
        file.write(json_study)
    with open("insert_into_data_files.json", "w") as file:
        file.write(json_data_files)

def main() -> None:
    """Module's entry-point."""

    fields = (
        'run_accession',
        'fastq_ftp',            # uses ftp, and contains both reads (if paired) separated by a semi colon (;)
        'fastq_md5',            # \ again, both separated by semicolon (;)
        'tax_id',
        'sample_accession',
        'study_accession',
        'library_strategy',     # read_type (it refers to wether a run is rna-seq, wgs, ...)
        'instrument_platform',  # platform
        'library_layout',       # paired
        'experiment_alias',     # \
        'experiment_title',     #  \ go to the same place
        'description',
        'library_source',       # library_name (name is requested later, but it's often empty so we'll catch it still to fill tissue-related info if available but best to use source for this field)
        'library_selection',    # is this really needed?
        'tissue_type',          # \
        'tissue_lib',           #  \ these both work towards filling the same field
        'cell_line',            # \ 
        'cell_type',            #  \ same as before
        'library_name',         # see library_source note
        'strain',
        'cultivar',             # \
        'ecotype',              #  \
        'isolate',              #   \ these last ones are to cater for the different ways of expresin "strains"
        # 'cage_protocol',        # this could be important to some farmed data, so added
        'center_name'
    )

    parser = argparse.ArgumentParser(prog='get_metadata.py',description='Get metadata for a given run accession from the ENA site')
    parser.add_argument(
        '--run',
        '-r',
        default="SRR4240445",
        help="Run accession to grab metadata from",
    )
    args = parser.parse_args()
    
    data = _request_data(args.run,fields)
    json = _json_parse(data,fields)

    #return json # not sure if this is it, or for nextflow to grab the output it has to be printed... if it needs printing, next line would solve it...
    print(json)

if __name__ == "__main__":
    main()
