


import requests
import os
import pymysql
import argparse

def parse_data(data):
    """Parse NCBI results and provides a useful dictionary containing all relevant information 

    Args:
        data (dict): json/dictionary obtained from NCBI API

    Returns:
        dict: dictionary with relevant data in a flatten format
    """
    
    ncbi_parsed = {
    'lowest_taxon_id' : data['reports'][0]['organism']['tax_id'],
    'biosample_id' : data['reports'][0]['assembly_info']['biosample']['accession'],
    'gca_chain' : data['reports'][0]['accession'].split('.')[0],
    'gca_version' : data['reports'][0]['accession'].split('.')[1],
    'is_current' : data['reports'][0]['assembly_info']['assembly_status'],
    'asm_type' : data['reports'][0]['assembly_info']['assembly_type'],
    'asm_level' : data['reports'][0]['assembly_info']['assembly_level'],
    'asm_name' : data['reports'][0]['assembly_info']['assembly_name'],
    'release_date' : data['reports'][0]['assembly_info']['release_date'],
    'submitter' : data['reports'][0]['assembly_info']['submitter'],
    'bioproject_id' : data['reports'][0]['assembly_info']['bioproject_accession'],
    'dtol_id' : '',
    'scientific_name' : data['reports'][0]['organism']['organism_name'],
    'dtol_prefix' : '',
    'clade' : '',
    'parlance_name' : '' ,
    'stable_id_prefix' : '',
    'assembly_metrics' : data['reports'][0]['assembly_stats']
    }
    
    try:
        common_name = data['reports'][0]['organism']['common_name']
    except:
        common_name = ""
    
    species_taxon_id = data['reports'][0]['organism']['tax_id']
    
    try:
        refseq_accession = data['reports'][0]['paired_accession']
    except:
        #print('No refseq accession found, set as NULL')
        refseq_accession = ''
        
    try:
        intra_name = list(data['reports'][0]['organism']['infraspecific_names'].keys())[0]
        intra_type = list(data['reports'][0]['organism']['infraspecific_names'].values())[0]  
    except:
        intra_name = ''
        intra_type = ''
    
    
    ncbi_parsed.update({'refseq_accession':refseq_accession,
                        'intra_type':intra_name,
                        'intra_name':intra_type,
                        'species_taxon_id':species_taxon_id,
                        'common_name' : common_name})
    
    return ncbi_parsed

def create_query(table_name,metadata, db_params, table_conf):
    """create mysql insert query for the specified table and retrieve the values from dictionary
    
    Args:
        table_name (str): name of the table from the assembly metadata database
        metadata (dict): dictionary with relevant data in a flatten format
        db_params (dict): host and credential to connect to db
        table_conf (dict): dictionary with the tables and their 'fill method' and primary key name
        
    Returns:
        str: mysql query to insert metadata in specified table
    """
    
    ## Connecting DB
    conn = pymysql.connect(
        host=db_params['host'],
        user=db_params['user'],
        password=db_params['password'],
        database=db_params['database'],
        port=int(db_params['port']))
    cur  = conn.cursor()
    
    ## Getting variables names: 
    cur.execute(f"DESCRIBE {table_name};")
    output = cur.fetchall()
    
    cur.close()
    conn.close()
    
    table_var = [] 
    for i in range(len(output)):
        # Exclude autoincremental key
        if output[i][5] != 'auto_increment': # and output[i][3] != 'MUL':
            table_var.append(output[i][0])
    table_var_string = ', '.join(table_var)
    
    method=table_conf[table_name]['method']
    
    # Getting value string
    if method=='per_col':
        #values_string = [metadata[var] for var in metadata if var in table_var]
        values_string = [metadata[x] for x in table_var]
        values_string = ", ".join([f"'{value}'" for value in values_string]).replace(", ''", ", NULL" )
        values_string = "({})".format(values_string)
        
    elif method =='per_row': #Method valid for assembly_metrics only
        value_list = []
        for key,value in metadata[table_name].items():
            value_item = f"('{metadata['assembly_id']}', '{key}', '{value}')"
            value_list.append(value_item)
            values_string =  ', '.join(value_list)
    else:
        print('Populate method not recognized')     
        

    # Writing query
    insert_query = f"""INSERT INTO {table_name} ({table_var_string}) VALUES {values_string}"""
    
    return insert_query

def execute_query(query, table_conf, table_name, metadata, db_params):
    """Execute insery query and 

    Args:
        query (str): mysql query to insert data
        table_conf (dict): dictionary with the tables and their 'fill method' and primary key name
        table_name (str): table name
        metadata (dict): dictionary with relevant data in a flatten format 
        db_params (dict): host and credential to connect to db

    Returns:
        dict: updated metadata for following queries that might need the auto_incremental/primary key 
        
    """
    
    conn = pymysql.connect(
        host=db_params['host'],
        user=db_params['user'],
        password=db_params['password'],
        database=db_params['database'],
        port=int(db_params['port']))
    cur  = conn.cursor()
    cur.execute(query)
    last_id = cur.lastrowid
    metadata.update({table_conf[table_name]['id_name'] : last_id})
    
    cur.close()
    conn.close()
    
    return metadata

def main() -> None:
    
    """Module's entry-point
    """
    
    db_params = {
        'host':'mysql-ens-genebuild-prod-1',
        'user':'ensadmin',
        'password':'ensembl',
        'port': '4527',
        'database' : 'gb_assembly_metadata_testing'
    }
    
    table_conf = {'assembly' : {'method': 'per_col',
                            'id_name': 'assembly_id'},
                'organism' : {'method': 'per_col',
                            'id_name': 'organism_ID'},
                'species'  : {'method': 'per_col',
                            'id_name': 'species_id'},
                'assembly_metrics' : {'method': 'per_row',
                                    'id_name': 'assembly_metrics_id'}
    }

    parser = argparse.ArgumentParser(prog='populate_metadata_tables.py', 
                                     description='Using the list of accession retrieves and save the metadata in the DB')
    
    parser.add_argument('--file-path', type=str,
                        help="File with the list of GCA accession to register in text plain format")
    
    args = parser.parse_args()
    print(args)
    
    if os.path.exists(args.file_path):

        with open(args.file_path, 'r') as file:
            for gca in file:
                print(gca.strip())
                
                uri = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/"+str(gca).strip()+'/dataset_report'
                response = requests.get(uri) #, params=ncbi_params)
                response.raise_for_status()
                data = response.json()
                
                metadata = parse_data(data)
                
                for table_name in table_conf:
                    #print(f"Populating {table_name} table")
                    query =  create_query(table_name, metadata, db_params, table_conf)
                    metadata = execute_query(query, table_conf, table_name, metadata, db_params)
        
        file.close()
    
    else:
        print("GCA list file not found. Re run to generate file ")







