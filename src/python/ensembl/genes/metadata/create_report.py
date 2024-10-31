import logging
import argparse
import pymysql
import json
import os

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
    with pymysql.connect(**metadata_params) as con:
        with con.cursor() as cur:
            cur.execute(query)
            return cur.fetchall()

def fetch_report_data(metadata_params, gca_list, bioprojects):
    gca_string = ','.join([f"'{gca}'" for gca in gca_list])
    bioprojects_string = ','.join([f"'{project}'" for project in bioprojects])

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
            SELECT assembly.GCA_chain, bioproject.bioproject_id, species.scientific_name, species.clade
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

def create_report(args, report_data):
    sections = [
        ("Num GCA", len(report_data['gca_list'])),
        ("Update Date", args.update_date),
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
        default=['PRJEB47820', 'PRJEB61747', 'PRJEB40665', 'PRJEB43743', 'PRJNA489243', 'PRJNA813333'],
        help='List of bioproject IDs to filter assemblies (default: specific bioprojects)'
    )

    args = parser.parse_args()
    logging.info(f"Script arguments: {args}")

    gca_list = load_file_lines(args.file_list)
    metadata_params = load_json(args.metadata)

    report_data = fetch_report_data(metadata_params, gca_list, args.bioprojects)
    report_data['gca_list'] = gca_list 

    create_report(args, report_data)
    logging.info("Report generated successfully.")

if __name__ == '__main__':
    main()