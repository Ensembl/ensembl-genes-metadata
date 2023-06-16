""" This script retrieves RefSeq assembly report summary from GenBank and compares the entries with the meta database.
    All assemblies in the meta database with an an updated RefSeq accession is updated """ 

import re
import urllib.request,json
from urllib.error import HTTPError
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import text
import sys
from Bio import Entrez
import requests
import os

Entrez.email = os.getenv('genebuild_email')

import pymysql
pymysql.install_as_MySQLdb()

def update_refseq_records(host,user,password,port,database,user_r):
    """ Retrieve RefSeq summary report and check against meta database entries """

    #Establish database connection
    con = pymysql.connect(
                host=host, user=user, passwd=password, port=port, database=database.strip(), charset='utf8mb4', cursorclass=pymysql.cursors.DictCursor
            ) 
    total_read = 0
    name_cnt = 0
    print('Downloading RefSeq records')
    ftp_path = 'http://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt'
    response = requests.get(ftp_path)
    header_pattern = '^(\w+.*)$'
    content_pattern = '^[GCA]'
    encoding_style = response.encoding
    assemblies_to_update = {}
    name_update = {}
    name_update_rpt = {}
    refseq_update_rpt = {}
    existing_accessions = {}
    acc_format = re.compile(r"^GC[AF]_\d+\.\d+$") #regex to filter genbank accession
            
    #Get all existing assembly ids
    sql = "SELECT CONCAT(chain,'.',version), assembly.assembly_id, refseq_accession, assembly_name FROM assembly, meta WHERE assembly.assembly_id = meta.assembly_id"
    results = fetch_data(sql,database,host,port,user_r,'',user_r)
    #Keep accessions and assembly ids in separate dictionaries
    for row in results:
        existing_accessions[row[0]] = row[1]
        name_update[row[0]] = row[3]
        if row[2] is None:
            assemblies_to_update[row[1]] = 1
        else:
            assemblies_to_update[row[1]] = row[2]

    try:
        for entry in urllib.request.urlopen(ftp_path):
            if (re.search(r'README_assembly_summary.txt',str(entry)) or re.search(r'assembly_accession',str(entry))):
                continue
            entry = entry.decode()
            line = entry.split('\t')
            line[0].strip
            # Check that there is a valid genbank accession for each refseq entry inorder to proceed with potential update run
            if (line[17]  != 'na'):
                line[17].strip
                if not re.match(acc_format, line[17]):
                    raise ValueError('Accession '+ line[17] + ' not in the right format ')
                if line[17] in existing_accessions:#RefSeq entry exists in database
                    if assemblies_to_update[existing_accessions[line[17]]] == 1:
                        #RefSeq accession requires updating
                        with con.cursor() as cur:
                            cur.execute('UPDATE meta SET refseq_accession = %s WHERE assembly_id = %s',  (line[0], existing_accessions[line[17]]))
                            #con.commit()
                        total_read += 1
                        #Store each update for reporting purposes
                        refseq_update_rpt[line[17]] = line[17] + '\t' + line[15] + '\t' + line[0] + '\n'
                    if name_update[line[17]] != line[15]:
                        #Assembly name requires updating
                        with con.cursor() as cur:
                            cur.execute('UPDATE meta SET assembly_name = %s WHERE assembly_id = %s',  (line[15], existing_accessions[line[17]]))
                            #con.commit()
                        name_cnt += 1
                        name_update_rpt[line[17]] = line[17] + '\t' + name_update[line[17]] + '\t' + line[15] + '\n'
    except urllib.error.HTTPError as e:
        if e.code == 404:
            print('No response was received. Possibly missing report file.\n')
    else:
        con.commit()
    slack_reporting(refseq_update_rpt,name_update_rpt)
    return

""" Report update status if any to the respective slack channel """
def slack_reporting(refseq_update_rpt,name_update_rpt):
    name_rpt = os.path.join(os.getenv('meta_database_dir'),'name_update.txt')
    rpt = os.path.join(os.getenv('meta_database_dir'),'refseq_summary_update.txt')
    refseq_report_line = "Genome accession\tAssembly name\tRefseq accession\n"
    name_report_line = "Genome accession\tOld assembly name\tNew assembly name\n"
    ref_cnt = 0
    asm_cnt = 0
    if name_update_rpt:
        for name_report in name_update_rpt:
            with open(name_rpt, 'a') as writer:
               writer.write(name_update_rpt[name_report])
               name_report_line = name_report_line + name_update_rpt[name_report]
            writer.close()
            asm_cnt += 1
    if refseq_update_rpt:
        for report in refseq_update_rpt:
            with open(rpt, 'a') as writer:
               writer.write(refseq_update_rpt[report])
               refseq_report_line = refseq_report_line + refseq_update_rpt[report]
            writer.close()
            ref_cnt += 1
    #Messages header
    name_report_line = 'Number of assemblies with assembly name update = ' + str(asm_cnt) + '\n' + name_report_line
    refseq_report_line = 'Number of assemblies with RefSeq accession update = ' + str(ref_cnt) + '\n' + refseq_report_line
    #Check what kind of update if any was done
    if ref_cnt > 0 and asm_cnt > 0:
        summary_report = refseq_report_line + name_report_line
    elif ref_cnt > 0 and asm_cnt == 0:
        summary_report = refseq_report_line
    elif ref_cnt == 0 and asm_cnt > 0:
        summary_report = name_report_line
    else:
        summary_report = 'No update to refseq or assembly name this time'

    #Set message token for reporting
    payload="{\"channel\": \"@denye\", \"username\": \"registry_messenger\", \"text\": \"" + summary_report  +"\"}"
    url = os.getenv('slack_token')
    headers = {'content-type': 'application/json', 'Accept-Charset': 'UTF-8'}
    r = requests.post(url, data=payload, headers=headers)

""" Method to access and query the meta database """
def fetch_data(query, database, host, port, user, password, user_r):
    """ This method fetches the requested data from the relevant database. """
    try:
        conn = pymysql.connect(
            host=host, user=user, passwd=password, port=port, database=database.strip()
        )

        cursor = conn.cursor()
        cursor.execute(query)
        info = cursor.fetchall()
    except pymysql.err.InternalError as e:
        print(str(e))
    except Exception as e:
        print(str(e))

    cursor.close()
    conn.close()
    return info  

""" Main method """
if __name__ == '__main__':
    """Main script entry-point."""
    update_refseq_records(os.getenv('GBS1'),os.getenv('GBUSER'),os.getenv('GBPASS'),int(os.getenv('GBP1')),os.getenv('REG_DB'),os.getenv('GBUSER_R'))
