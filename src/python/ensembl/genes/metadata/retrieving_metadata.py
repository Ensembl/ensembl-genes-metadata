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

"""This module retrieves metadata from the NCBI API for a given GCA accession 
        and stores it in JSON files to be inserted in the database.
        
    Args:
        accession (str): GCA accession to retrieve metadata

    Returns:
        str: json file with assembly's metadata
        str: json-like file with metrics and bioproject lineage with .tmp extension
        str: json-like file with species data with .tmp extension
"""

import requests
import argparse
import json
import logging
from tenacity import retry, stop_after_attempt, wait_random
from typing import Dict, List

@retry(stop=stop_after_attempt(5), wait=wait_random(min=1, max=10))
def connection_ncbi(uri: str) -> requests.Response:
    response = requests.get(uri)
    response.raise_for_status()
    return response

def parse_data(data: Dict) -> List[Dict]:
    """Parse NCBI results and provides a useful dictionary containing all relevant information 
    
    Args:
        data (dict): original report obtained from NCBI API
    
    Returns:
        list[dict]: a list of three dictionaries with relevant metadata to store in db
    """
    
    data_dict = {'assembly': {}}
    
    data_dict_2 = {
        'organism': {},
        'assembly_metrics': {},
        'bioproject': {}
    }
    
    data_dict_3 = {
        'species': {}
    }
    
    # Setting Optional keys first
    try:
        common_name = data['reports'][0]['organism']['common_name'].replace("'", "''")
    except:
        common_name = ""
    
    try:
        refseq_accession = data['reports'][0]['paired_accession']
    except:
        refseq_accession = ""
        
    try:
        infra_type = list(data['reports'][0]['organism']['infraspecific_names'].keys())[0]
        infra_name = list(data['reports'][0]['organism']['infraspecific_names'].values())[0].replace("'", "''")  
        if infra_type == 'sex':
            logging.info("Assembly have incorrect infraspecific type: sex. Setting to empty")
            infra_type = ""
            infra_name = ""
    except:
        infra_name = ""
        infra_type = ""
        
    try:
        biosample_id = data['reports'][0]['assembly_info']['biosample']['accession']
    except:
        biosample_id = ""
    
    # Building dictionaries for assembly metadata tables
    data_dict['assembly'].update({
        'lowest_taxon_id' : data['reports'][0]['organism']['tax_id'],
        'gca_chain' : data['reports'][0]['accession'].split('.')[0],
        'gca_version' : data['reports'][0]['accession'].split('.')[1],
        'is_current' : data['reports'][0]['assembly_info']['assembly_status'],
        'asm_type' : data['reports'][0]['assembly_info']['assembly_type'],
        'asm_level' : data['reports'][0]['assembly_info']['assembly_level'],
        'asm_name' : data['reports'][0]['assembly_info']['assembly_name'],
        'refseq_accession': refseq_accession,
        'release_date' : data['reports'][0]['assembly_info']['release_date'],
        'submitter' : data['reports'][0]['assembly_info']['submitter'].replace("'", "''").lstrip('\ufeff')
    })
    
    data_dict_3['species'].update({
        'lowest_taxon_id' : data['reports'][0]['organism']['tax_id'],
        'scientific_name' : data['reports'][0]['organism']['organism_name'].replace("'", "''") ,
        'common_name' : common_name
    })
    
    data_dict_2['organism'].update({
        'biosample_id' : biosample_id,
        'bioproject_id' : data['reports'][0]['assembly_info']['bioproject_accession'],
        'infra_type':infra_type,
        'infra_name':infra_name
    })
    
    data_dict_2.update({'assembly_metrics':data['reports'][0]['assembly_stats']})
    
    # Parsing bioproject metadata 
    bioproject_lineage = {} 
    seen_accessions = set()
    
    bioproject_dict = data['reports'][0]['assembly_info']['bioproject_lineage'][0]['bioprojects']
    for item in bioproject_dict:
        accession = item['accession']
        title = item['title'].replace("'", "")
        # Check if accession is not seen before
        if accession not in seen_accessions:
            seen_accessions.add(accession)
            bioproject_lineage[accession] = title
    
    data_dict_2.update({'bioproject':bioproject_lineage})
    
    return data_dict, data_dict_2, data_dict_3

def main():
    """Module's entry-point
    """
    logging.basicConfig(filename="retrieving_metadata.log", level=logging.DEBUG, filemode='w',
                    format="%(asctime)s:%(levelname)s:%(message)s")
    
    parser = argparse.ArgumentParser(prog='retrieving_metadata.py', 
                                    description="Retrieve metadata from NCBI API for a given GCA accession and store it in JSON files to be inserted in the database.")
    
    parser.add_argument('--accession', 
                        type=str,
                        help='GCA accession to retrieve metadata')
    
    args = parser.parse_args()
    logging.info(args)
    accession = args.accession.strip()
    
    uri = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{accession}/dataset_report"
    response = connection_ncbi(uri)
    data = response.json()
    
    logging.info(f"Retrieved data for {accession}")
    dict1, dict2, dict3 = parse_data(data)
    
    logging.info("Saving data in JSON files")
    file1 = f"{accession}_assembly.json"
    with open(file1, 'w') as file:
        json.dump(dict1, file)
    file.close()
    
    file2 = f"{accession}_metadata.tmp"
    with open(file2, 'w') as file:
        json.dump(dict2, file)
    file.close()
    
    file3 = f"{accession}_species.tmp"
    with open(file3, 'w') as file:
        json.dump(dict3, file)
    file.close()

if __name__ == '__main__':
    main()