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
"""
Get metadata for transcriptomic registry module

This script simply gathers all metadata required to be added to out 
transcriptomic database, parses it from the ENA database, and returns a json
to be parsed into the final DB.

parameters
----------
arg1    : string
    string with run_accession to search

returns
-------
stout   : string
    contains a simplified json in the agreed upon format to fill the DB
"""

import argparse
import requests

def _request_data(run_accession) -> str:
    """
    Request the data from ENA site

    parameters
    ----------
    accession_run : string
        accession id for run

    returns
    -------
    tsv_list : string
        list of requeste data for a given accession_run
    """

    query = 'run_accession="' + run_accession + '"'

    print(query)

    data = {
    'result': 'read_run',
    'query': query,
    'fields': 'run_accession,center_name,sample_accession,description,cell_line,cell_type,cultivar,ecotype,experiment_alias,experiment_title,isolate,library_name,sample_alias,sample_description,serotype,strain,study_accession,tax_id,tissue_lib,tissue_type',
    'format': 'tsv',
    }

    response = requests.post('https://www.ebi.ac.uk/ena/portal/api/search', data=data)
    return response

#def _json_parse()

def main() -> None:
    """Module's entry-point."""
    print("Parsing arguments...")
    parser = argparse.ArgumentParser(prog='get_metadata.py',description='Get metadata for a given run accession from the ENA site')
    parser.add_argument(
        "run_accession",
        type=str,
        help="Run accession to grab metadata from",
    )
    args = parser.parse_args()
    print("Requesting data...")
    data = _request_data(args.run_accession)
#    json = _json_parse(data)
