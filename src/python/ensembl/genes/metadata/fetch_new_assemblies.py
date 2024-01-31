


"""Module will connect to NCBI API and Genebuild databases to determinate the list of 
GCA accession to register in the database

Returns:
    file: txt file containing a list of accession
"""

import requests
import os
from datetime import datetime
import pymysql
import argparse

def set_date(taxon, ncbi_params, file_path):
    """Set a date to retrival new assemblies from NCBI

    Args:
        file_path (str): path with the last database update (optional) 
        taxon (int): Taxon ID, by default it is 3745 (Eukaryote domain)
        ncbi_params (dict): NCBI API's parameters

    Returns:
        dict: Updated NCBI API's parameters
        str: date to be used to retrieve assemblies 
    """
    print('Checking date to retrive assemblies')
    
    release_date = "01/01/2019"

    if os.path.exists(file_path):
        try:
            with open(file_path, 'r') as file:
                for line in file:
                    if int(line.split()[0]) == taxon:
                        release_date = line.split()[1]
                        ncbi_params['filters.first_release_date']= release_date
                        print(f"Last update record found! {line.split()[1]} will be used to retrive new available assemblies for taxon {taxon}")
                        break
                else:
                    print(f"There is not last update record saved for taxon {taxon}. Using {release_date} as default date to retrive assemblies")

        except (ValueError, IndexError):
            print (f"Error: Unable to parse the file. Using {release_date} as default date to retrive assemblies")
    else:
        print(f"Last update file not found. Using {release_date} as default date to retrive assemblies")
  
    return ncbi_params, release_date

def fetch_gca_list (taxon, ncbi_params):
    """Fetch a list of accessions avilable since the last update date 

    Args:
        taxon (int): Taxon ID, by default it is 3745 (Eukaryote domain)
        ncbi_params (dict): NCBI API's parameters

    Returns:
        set: a list of GCA accession
    """
    gca_list = [] 
    page_token=None
    
    ## build URI
    try:
        int(taxon)
        uri = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon/"+str(taxon)+'/dataset_report'

    except ValueError:
        print("Error: Please enter a valid taxon numeric value.")
        
    except TypeError:
        print(f"Error: Invalid taxon value to build query: {taxon}")
        
    else:
        print(f"Valid value to connect to NCBI API and access to genome reports per taxon: {taxon}")

    # Connecting to API and Fetch GCA list
    while True:

        if page_token is not None: 
            #print('Reading Next page!')
            ncbi_params['page_token']=response.headers['X-Ncbi-Next-Page-Token']
        
        try:
            response = requests.get(uri, params=ncbi_params) # Sending query
            # print(response.url)
            response.raise_for_status()  # Raises an HTTPError if the HTTP request returned an unsuccessful status code
            #print(response.url)
            data = response.json()
            
            N = len(data['reports'])
            gca_list.extend([data['reports'][i]['accession'] for i in range(N)])
            
            page_token=response.headers['X-Ncbi-Next-Page-Token']
    
        # Error request
        except requests.exceptions.RequestException as e:
            print(f"Error Request: {e}")
        
        # Error in reading token, there are not more pages.
        except KeyError:
            print('All available assemblies were retrived')
            break

    return set(gca_list)       

def build_db_query(release_date):
    """ Build mysql queries to retrieve data from the databases
    
    Args:
        release_date (str): date to be used to retrieve assemblies
    
    Returns:
        dict: a dictionarya containing a query for each database
    
    """
    
    release_date_sql = datetime.strptime(release_date, '%m/%d/%Y').strftime('%Y-%m-%d')
    
    query_registry =  """SELECT concat(assembly.chain, '.', assembly.version) AS GCA
    FROM assembly
    INNER JOIN meta
    ON assembly.assembly_id = meta.assembly_id
    WHERE meta.assembly_date >= DATE('{}') ; """.format(release_date_sql)
    
    query_metadata = """SELECT concat(gca_chain, '.', gca_version) 
    FROM assembly 
    WHERE release_date >= DATE('{}') ; """.format(release_date_sql)
    
    db_query = {
    'asm_metadata':
        {'db': 'gb_assembly_metadata_testing',
        'query': query_metadata},
    'asm_registry':
        {'db': 'gb_assembly_registry',
        'query': query_registry},
    }

    return db_query

def fetch_records_db(db_params, database, query):
    """Fetch assemblies that have been register after the last update

    Args:
        db_params (dict): NCBI API's parameters
        database (str): database to be used to fetch records
        query (_type_): query for the database 

    Returns:
        list: list of GCA recorded after the last update date
    """
    
    try:
    
        ## Connection to the DB
        
        conn = pymysql.connect(
            host=db_params['host'],
            user=db_params['user'],
            database=database,
            port=int(db_params['port']))
        cur  = conn.cursor()
        # Querying database
        cur.execute(query)
        output= cur.fetchall() 
        
        # Processing results
        reg_gca = [output[i][0] for i in range(len(output))]
    
    except pymysql.err.InternalError as e:
        print(str(e))
        
    except Exception as e:
        print(str(e))  
        
    finally:
        conn.close()
        
    return reg_gca

def get_gca_register(db, db_query, gca_list, output_path, db_params):
    """ Get list of GCA assemblies to register

    Args:
        db (str): option indicating what database use to retrive accession and later compare
                    (asm_metadata, asm_registry, both)
        db_query (dict): database and its query
        gca_list (list): list of GCA accessions retrived from API
        output_path (str): Path to save output
    """
    
    output_file = str(output_path) + "assemblies_to_register.txt"

    
    accessions_records = []

    if db in ['asm_metadata', 'asm_registry']:
        print(f'Assemblies will be updated according to {db} records only')
        accessions_records = fetch_records_db(db_params=db_params, database=db_query[db]['db'], query=db_query[db]['query'] ) 
        
        
    elif db == 'both':

        print('Assemblies will be updated using the records from both databases')
        for db in db_query:
            accessions_records = accessions_records + fetch_records_db(db_params=db_params, database=db_query[db]['db'], query=db_query[db]['query'] ) 
            
    accessions_to_register = gca_list - set(accessions_records)
    
    with open(output_file, 'w') as file:
        file.write('\n'.join(accessions_to_register))
    file.close()
        
    print(f'Accessions to register: {len(accessions_to_register)}')
    
def main():
    """module's entry-point
    """
    ncbi_params = {
                'bioprojects': None,
                'filters.reference_only': 'false',
                'filters.assembly_source': 'genbank',
                'filters.has_annotation': 'false',
                'filters.exclude_paired_reports': 'false',
                'filters.exclude_atypical': 'true',
                'filters.assembly_version': 'all_assemblies',
                'filters.assembly_level': ['scaffold', 'chromosome', 'complete_genome', 'contig'],
                'filters.first_release_date': '01/01/2019',
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

    db_params = {
        'host':'mysql-ens-genebuild-prod-1',
        'user':'ensro',
        'port': '4527'
    }
            
    parser = argparse.ArgumentParser(prog='fetch_new_assemblies.py', 
                                     description='Determinate which assemblies register')
    
    parser.add_argument('--taxon', 
                        default=3745,
                        type=int,
                        help='Valid Taxon id: eukaryota - 3745')
    parser.add_argument('--file-path',
                        default='',
                        type=str,
                        help="File with the last update date", 
                        action='append')
    parser.add_argument('--output-path',
                        type=str,
                        help="Output path, file will be named as assemblies_to_register.txt", 
                        action='append')
    
    args = parser.parse_args()
    print(args)
    
    ncbi_params, release_date = set_date(args.taxon, ncbi_params, args.file_path)
    gca_list = fetch_gca_list(args.taxon, ncbi_params)
    db_query = build_db_query(release_date)
    get_gca_register('both', db_query, gca_list, args.output_path, db_params)

if __name__ == '__main__':
     main()