## Pending: 
## - bioproject
## - retry connection to NCBI API
## - choose a better raise exception 

"""_summary_

    Raises:
        ValueError: _description_

    Returns:
        _type_: _description_
"""

import requests
import argparse
import json

def parse_data(data):
    """Parse NCBI results and provides a useful dictionary containing all relevant information 
    
    Args:
        data (dict): json/dictionary obtained from NCBI API
    
    Returns:
        list[dict]: A list of two dictionaries with relevant metadata to store in db
    """
    
    data_dict = {
    'assembly': {},
    'organism': {}
    }
    
    data_dict_2 = {
        'assembly_metrics': {}
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
        infra_name = list(data['reports'][0]['organism']['infraspecific_names'].keys())[0].replace("'", "''")
        infra_type = list(data['reports'][0]['organism']['infraspecific_names'].values())[0]  
    except:
        infra_name = ""
        infra_type = ""
    
    # Building dictionaries for assembly metadata tables
    data_dict['assembly'].update({
        'lowest_taxon_id' : data['reports'][0]['organism']['tax_id'],
        'biosample_id' : data['reports'][0]['assembly_info']['biosample']['accession'],
        'gca_chain' : data['reports'][0]['accession'].split('.')[0],
        'gca_version' : data['reports'][0]['accession'].split('.')[1],
        'is_current' : data['reports'][0]['assembly_info']['assembly_status'],
        'asm_type' : data['reports'][0]['assembly_info']['assembly_type'],
        'asm_level' : data['reports'][0]['assembly_info']['assembly_level'],
        'asm_name' : data['reports'][0]['assembly_info']['assembly_name'],
        'refseq_accession': refseq_accession,
        'release_date' : data['reports'][0]['assembly_info']['release_date'],
        'submitter' : data['reports'][0]['assembly_info']['submitter']
    })
    
    data_dict_3['species'].update({
        'lowest_taxon_id' : data['reports'][0]['organism']['tax_id'],
        'scientific_name' : data['reports'][0]['organism']['organism_name'],
        'common_name' : common_name
    })
    
    data_dict['organism'].update({
        'biosample_id' : data['reports'][0]['assembly_info']['biosample']['accession'],
        'bioproject_id' : data['reports'][0]['assembly_info']['bioproject_accession'],
        'infra_type':infra_name,
        'infra_name':infra_type
    })
    
    data_dict_2.update({'assembly_metrics':data['reports'][0]['assembly_stats']})
    
    # Parsing bioproject metadata 
    bioproject_lineage = {} ### Create empty dictionary to later add to 
    seen_accessions = set() ### Need to 
    
    bioproject_dict = data['reports'][0]['assembly_info']['bioproject_lineage'][0]['bioprojects']
    for item in bioproject_dict:
        accession = item['accession'] # bioproject accession PRJEBXXXX
        title = item['title'].replace("'", "")
        # Check if accession is not seen before
        if accession not in seen_accessions:
            seen_accessions.add(accession)
            bioproject_lineage[accession] = title
    
    #data_dict_2.update({'bioproject_lineage':bioproject_lineage})
    
    return data_dict, data_dict_2, data_dict_3

def main():
    
    parser = argparse.ArgumentParser(prog='retrieving_metadata.py', 
                                    description="Retrieves accession's metadata and store it in JSON files")
    
    parser.add_argument('--accession', 
                        default="GCA_036172605.1",
                        type=str,
                        help='Valid GCA accession')
    
    parser.add_argument('--output-path',
                        default="/Users/vianey/Documents",
                        type=str,
                        help="Output path, file will be named with the accession as prefix")
    
    args = parser.parse_args()
    accession = args.accession.strip()
    
    attempt_num = 1
    success=True
    """
    while success == True or attempt_num<=3:
        uri = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{str(args.accession).strip()}/dataset_report"
        response = requests.get(uri)
        response.raise_for_status()
        
        if response.raise_for_status()is None:
            success = True
        else:
            attempt_num+=1
    """
    uri = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{accession}/dataset_report"
    response = requests.get(uri)
    response.raise_for_status()
    
    if success:
        data = response.json()
        
        dict1, dict2, dict3 = parse_data(data)
        
        # Save first json for insert
        #file1 = f"{args.output_path}/{args.accession}_metadata.json"
        file1 = f"{accession}_metadata.json"
        with open(file1, 'w') as file:
            json.dump(dict1, file)
        file.close()
        
        # Save second json for next module
        #file2 = f"{args.output_path}/{args.accession}_metrics.json"
        file2 = f"{accession}_metrics.tmp"
        with open(file2, 'w') as file:
            json.dump(dict2, file)
        file.close()
        
        file3 = f"{accession}_species.tmp"
        with open(file3, 'w') as file:
            json.dump(dict3, file)
        file.close()
    else:
        raise ValueError(f"It was not possible to connect to NCBI API datasets for {accession}")

if __name__ == '__main__':
    main()