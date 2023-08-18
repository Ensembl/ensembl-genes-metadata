# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

""" Script to check for the availability of transcriptomic data when provided with the taxon_id """

import sys
import argparse
import re
import requests
import pymysql
pymysql.install_as_MySQLdb()

#Method to establish connection to database and query it
def fetch_data(query, database, host, port, user, password):
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

""" Use taxon_id to check for the existence of short reads (paired ended) or long reads from the ENA """
def transcriptomic_status(output_path,host,user,password,port,db):
    """ Check status of transcriptomic data both in the registry for exisiting species and publicly for new species """
    ena_base_url = 'https://www.ebi.ac.uk/ena/portal/api/search?display=report'
    files_domain = 'domain=read&result=read_run'
    files_fields = 'run_accession,study_accession,experiment_accession,sample_accession,secondary_sample_accession,instrument_platform,instrument_model,library_layout,library_strategy,nominal_length,read_count,base_count,fastq_ftp,fastq_aspera,submitted_ftp,fastq_md5,library_source,library_selection,center_name,study_alias,experiment_alias,experiment_title,study_title'
    download_method = 'ftp'
    total_read = 0
    entry = 0
    status = 'not available'

    #Set db connection settings
    con = pymysql.connect(
            host=host, user=user, passwd=password, port=int(port), db=db.strip(), charset='utf8mb4', cursorclass=pymysql.cursors.DictCursor
    )
    #get species with rnaseq data set as not available
    sql = "SELECT distinct(species_id) FROM assembly WHERE rnaseq_data = 'not available'"
    info = fetch_data(sql,db,host,int(port),user,password)
    results = 0
    #Process each species for new transcriptomic data
    for taxon_id in info:
        print('Processing taxonomy '+str(taxon_id[0]))
        ftp_path = ena_base_url + '&query="tax_tree(' + str(taxon_id[0]) + ') AND instrument_platform=ILLUMINA AND library_source=TRANSCRIPTOMIC"&' + files_domain + '&fields=' + files_fields
        response = requests.get(ftp_path)
        fastq_file = 'fastq_' + download_method
        header_pattern = '^(\w+.*)$'
        content_pattern = '^[a-z]'
        encoding_style = response.encoding
        unwanted_reads = []
        wanted_reads = []

        if response:
            read_file = output_path + "/ena_read_set.txt"
            open(read_file, "wb").write(response.content)
            with open(read_file) as f:
                for line in f.readlines():
                    header = re.search(header_pattern, line)
                    if header:
                        read_entry = header.group(1)
                    content = re.search(content_pattern, read_entry)
                    if not content:
                        total_read += 1
                        wanted_reads.append(read_entry)
                        entry += 1
            if entry > 0:
                status = 'available'
                try:
                    with con.cursor() as cur:
                        cur.execute('UPDATE assembly SET rnaseq_data = %s WHERE species_id = %s' ,  (status, taxon_id))
                except Exception as error:
                    print ('Oops! Rnaseq status could not be updated. Kindly check the error below\n'+str(error))
                    break
                else:
                    results += 1
                    con.commit()
    return 'Number of species without data = '+str(len(info))+' and number found with data = '+str(results)

""" Main method to handle arguments """
if __name__ == '__main__':
    """Main script entry-point."""
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_path', help='Path to store transcriptomic data report', required=True)
    parser.add_argument('--port', help='Port number for host', required=True)
    parser.add_argument('--password', help='Password for host', required=False)
    parser.add_argument('--user', help='Mysql user', required=True)
    parser.add_argument('--dbname', help='Database to be backed up', required=True)
    parser.add_argument('--server', help='Host server for database', required=True)
    args = parser.parse_args()
    output_path = args.output_path
    host = args.server
    db = args.dbname
    user = args.user
    port = args.port
    password = args.password

    status = transcriptomic_status(output_path,host,user,password,port,db)
    print(status)
