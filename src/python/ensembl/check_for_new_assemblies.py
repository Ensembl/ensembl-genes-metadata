""" This script searches the public archives (NCBI) for all available eukaryotes and returns a list of unregistered assemblies """

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

def fetch_genbank_records():
    """ This method fetches publicly avaliable assembly records all all Eukaroytes. The records fetched and then formatted to be used by the registration script. It also returns the number as assemblies valid for registration """

    total_read = 0
    print('Downloading genbank records from 2017 onwards')
    ftp_path = 'http://ftp.ncbi.nlm.nih.gov//genomes/GENOME_REPORTS/eukaryotes.txt'
    response = requests.get(ftp_path)
    header_pattern = '^(\w+.*)$'
    content_pattern = '^[GCA]'
    encoding_style = response.encoding
    accessions_to_register = {}
    rpt = os.path.join(os.getenv('HOME'),'assemblies_to_register.txt')
    existing_accessions = {}
    acc_format = re.compile(r"^GC[AF]_\d+\.\d+$") #regex to filter genbank accession
            
    #Get all existing assembly ids
    sql = "SELECT CONCAT(chain,'.',version), assembly_id FROM assembly"
    results = fetch_data(sql,os.getenv('REG_DB'),os.getenv('GBS1'),int(os.getenv('GBP1')),os.getenv('GBUSER_R'),'')
    for row in results:
        existing_accessions[row[0]] = row[1]
    date_val = 0
    #Formatting output ready for use later in the registration module
    accessions = '['
    comma_replacement = ']'
    try:
        for entry in urllib.request.urlopen(ftp_path):
            if (re.search(r'#Organism/Name',str(entry))):
                continue
            entry = entry.decode()
            line = entry.split('\t')
            line[8].strip
            # We set a limit to retrieve data from 2017 upwards partly due to timing around which the major projects we collaborate with started. Also, most assemblies before 2017 were not of much higher quality either.
            if (int(line[14][:4]) >= 2017):
                date_val += 1
                if not re.match(acc_format, line[8]):
                    raise ValueError('Accession '+ line[8] + ' not in the right format ')
                if line[8] not in existing_accessions:
                    accessions_to_register[line[8]] = line[1]
                    total_read += 1
                    accessions = accessions + line[8] + ","
        accessions = comma_replacement.join(accessions.rsplit(',',1))
        with open(rpt, 'a') as writer:
           writer.write(accessions)
        writer.close()
    except urllib.error.HTTPError as e:
        if e.code == 404:
            print('No response was received. Possibly missing report file.\n')
    return total_read


def fetch_data(query, database, host, port, user, password):
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

if __name__ == '__main__':
    """Main script entry-point."""
    fetch_genbank_records()
