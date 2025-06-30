import logging
import argparse
import pymysql
import json
import os
from datetime import date 

def load_json(filepath):
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"{filepath} does not exist")
    with open(filepath, 'r') as f:
        return json.load(f)

def load_file_lines(filepath):
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"{filepath} does not exist")
    with open(filepath, 'r') as file:
        return [line.strip() for line in file if line.strip()]

def execute_query(metadata_params, query):
    logging.info(f"QUERY: {query}")
    conn = pymysql.connect(**metadata_params)
    with conn.cursor() as cur:
        cur.execute(query)
        return cur.fetchall()

def fetch_report_data(metadata_params, gca_list, bioprojects):
    gca_string = ','.join([f"'{gca}'" for gca in gca_list])
    bioprojects_string = ','.join([f"'{project}'" for project in bioprojects])
    logging.info(f"GCA string for report:{gca_string}")
    

    queries = {
        'asm_type': f"""
            SELECT asm_type, COUNT(*)
            FROM assembly
            WHERE CONCAT(GCA_chain, '.', gca_version) IN ({gca_string})
            GROUP BY asm_type;
        """,
        'asm_level': f"""
            SELECT asm_level, COUNT(*)
            FROM assembly
            WHERE CONCAT(GCA_chain, '.', gca_version) IN ({gca_string})
            GROUP BY asm_level;
        """,
        'bioproject_species': f"""
            SELECT CONCAT(assembly.GCA_chain, '.', assembly.gca_version), bioproject.bioproject_id, species.scientific_name
            FROM assembly
            INNER JOIN bioproject ON assembly.assembly_id = bioproject.assembly_id
            INNER JOIN species ON assembly.lowest_taxon_id = species.lowest_taxon_id
            WHERE bioproject.bioproject_id IN ({bioprojects_string})
            AND CONCAT(GCA_chain, '.', gca_version) IN ({gca_string});
        """,
        'refseq': f"""
            SELECT CONCAT(GCA_chain, '.', gca_version) AS GCA
            FROM assembly
            WHERE CONCAT(GCA_chain, '.', gca_version) IN ({gca_string})
            AND refseq_accession IS NOT NULL;
        """,
        'invalid_taxon': """
            SELECT CONCAT(GCA_chain, '.', gca_version) AS GCA, species.scientific_name, species.lowest_taxon_id, assembly.release_date
            FROM assembly
            INNER JOIN species ON assembly.lowest_taxon_id = species.lowest_taxon_id
            WHERE species_taxon_id = 0;
        """,
        'missing_biosample': """
            SELECT CONCAT(GCA_chain, '.', gca_version) AS GCA
            FROM assembly
            INNER JOIN organism ON assembly.assembly_id = organism.assembly_id
            WHERE organism.biosample_id IS NULL;
        """,
        'missing_submitter': """
            SELECT CONCAT(GCA_chain, '.', gca_version) AS GCA
            FROM assembly
            WHERE submitter IS NULL;
        """
    }

    report_data = {}
    report_data['gca_list'] = gca_list
    # Populate the rest of the report
    for key, query in queries.items():
        report_data[key] = execute_query(metadata_params, query)
    return report_data

def format_section(header, data):
    # Check if data is iterable (like list, tuple, etc.); if not, convert it to string directly
    if isinstance(data, (list, tuple)):
        formatted_data = "\n".join(str(row) for row in data)
    else:
        formatted_data = str(data)
    
    return f"{header}:\n{formatted_data if formatted_data else 'No data available'}\n\n"

def create_report(update_date, report_data):
    sections = [
        ("Num GCA", len(report_data['gca_list'])),
        ("Update Date", update_date),
        ("Overview Assembly Type", report_data['asm_type']),
        ("Overview Assembly Level", report_data['asm_level']),
        ("Assemblies from Relevant Bioprojects", report_data['bioproject_species']),
        ("Assemblies with RefSeq Available", report_data['refseq']),
        ("Assemblies wit Invalid taxons ID", report_data['invalid_taxon']),
        ("Missing Biosample", report_data['missing_biosample']),
        ("Missing Submitter", report_data['missing_submitter'])
    ]

    report = "\n".join([format_section(header, data) for header, data in sections])

    with open("report.txt", 'w') as file:
        file.write(report)

def create_csv_busco(gca_list, bioprojects, metadata_params):
    gca_string = ','.join([f"'{gca}'" for gca in gca_list])
    bioprojects_string = ','.join([f"'{project}'" for project in bioprojects])
        
    get_data_query = f"""
    select DISTINCT CONCAT(assembly.gca_chain, '.', assembly.gca_version) , assembly.lowest_taxon_id from assembly
    JOIN bioproject on assembly.assembly_id = bioproject.assembly_id
    JOIN taxonomy on taxonomy.lowest_taxon_id = assembly.lowest_taxon_id
    WHERE ((taxonomy.taxon_class_id = '40674' and assembly.lowest_taxon_id != '9606') 
    OR bioproject_id IN ({bioprojects_string}))
    and assembly.lowest_taxon_id != '9606' 
    and asm_type = 'haploid' and assembly.is_current = 'current' 
    and CONCAT(assembly.GCA_chain, '.', assembly.gca_version) IN ({gca_string});
    """
    output = execute_query(metadata_params, get_data_query)
    
    with open("gca_to_run_ncbi.txt", 'w') as file:
        file.write("gca,taxon_id\n")
        for gca, taxon_id in output:
            file.write(f"{gca},{taxon_id}\n")
        
def update_date(metadata_params):
    current_date = date.today()
    logging.info(f"Store last update date ({current_date}) in assembly metadata DB")
    try:
        connection = pymysql.connect(**metadata_params)
        with connection:
            with connection.cursor() as cursor:
                update_query = """UPDATE update_date SET date_value = %s WHERE update_type = 'regular_update';"""
                cursor.execute(update_query, (current_date,))
                logging.info("Update Successful")
    except Exception as e:
        logging.error(f"Failed to update date: {e}")
        raise  

def main():
    logging.basicConfig(
        filename="create_report.log", 
        level=logging.DEBUG, 
        filemode='w',
        format="%(asctime)s:%(levelname)s:%(message)s"
    )

    parser = argparse.ArgumentParser(
        prog='create_report.py',
        description='Create a report of the recently registered assemblies'
    )
    parser.add_argument('--file-list', help='Path to txt file with a list of GCA accessions', required=True)
    parser.add_argument('--metadata', help='Path to metadata database connection parameters', required=True)
    parser.add_argument('--update-date', help='Date used to update database', required=True)
    parser.add_argument(
        '--bioprojects', 
        nargs='+', 
        default=['PRJNA730823','PRJEB40665', 'PRJEB43510', 'PRJEB43743', 'PRJEB47820', 'PRJEB61747', '','PRJNA489243', 'PRJNA533106','PRJNA813333'],
        help='List of bioproject IDs to filter assemblies (default: specific bioprojects)'
    )

    args = parser.parse_args()
    logging.info(f"Script arguments: {args}")

    gca_list = load_file_lines(args.file_list)
    metadata_params = load_json(args.metadata)

    report_data = fetch_report_data(metadata_params, gca_list, args.bioprojects)
    report_data['gca_list'] = gca_list 

    create_report(args.update_date, report_data)
    logging.info("Report generated successfully.")
    
    logging.info("Create CSV file for Busco genome ")
    create_csv_busco(gca_list, args.bioprojects, metadata_params)

    update_date(metadata_params)

if __name__ == '__main__':
    main()