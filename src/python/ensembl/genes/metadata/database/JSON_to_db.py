import json
import pymysql
import argparse

def create_query(db_params, data, table_name='file_metrics'):
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
    
    ## Getting values to populate table
    value_list = []
    for accession in data.items():
        #print(accession[0])
        for metrics, status in accession[1].items():
            value_item = f"('{accession[0]}', '{metrics}', '{status['status']}')"
            value_list.append(value_item)
            values_string =  ', '.join(value_list)
    
    # Create Query       
    query = f"""INSERT INTO {table_name} ({table_var_string}) VALUES {values_string}"""
    
    return query

def execute_query(query, db_params):
    conn = pymysql.connect(
        host=db_params['host'],
        user=db_params['user'],
        password=db_params['password'],
        database=db_params['database'],
        port=int(db_params['port']))
    cur  = conn.cursor()

    try:
        cur.execute(query)
    except Exception as e:  
        print(f'Error: {e}')

    cur.close()
    
def main():

    db_params = {
        'host':'mysql-ens-genebuild-prod-1',
        'user':'ensadmin',
        'password':'ensembl',
        'port': '4527',
        'database' : 'transcriptomic_hackathon_vianey'
    }
    
    parser = argparse.ArgumentParser(prog='JSON_to_db.py', 
                                     description='Write FASTQC results to DB')
    
    parser.add_argument('--file-path', type=str,
                        help="JSON file containing fastqc results")
    
    args = parser.parse_args()
    
    file = open(str(args.file_path))
    data = json.load(file)
    file.close()
    
    query = create_query(db_params, data)
    execute_query(query, db_params)
    
if __name__ == '__main__':
     main()