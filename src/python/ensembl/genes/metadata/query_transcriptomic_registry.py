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

""" Script to retrieve candidate assemblies for annotation from the registry.  """

import pymysql
pymysql.install_as_MySQLdb()
import urllib.request,json
import argparse
import time
import sqlalchemy as db
import string
import re
import errno
import os
from urllib.error import HTTPError
from Bio import Entrez
Entrez.email = os.getenv('genebuild_email')

""" Check what task needs to be done and return the results """
def generate_report(db,host,port,user,output_file,taxon_id,task,task_list):
    results = 0
    if task in task_list.keys():
        print('Requested task ('+task+') in progress')
    else:
        print('Requested task cannot be performed. Check that it is included in the allowed list below')
        for tasks in task_list.keys():
            print('  '+tasks+' => '+task_list[tasks])
        quit()
    # Execute task specific queries and write results to out put file
    if task == 'get_short_read_at_species_level':
        sql = "SELECT * FROM short_read_data WHERE species_id = "+taxon_id
        info = fetch_data(sql,db,host,int(port),user,'')
        for row in info:
            results += 1
            filename = row[3].split(';')
            url = row[6].split(';')
            md5 = row[7].split(';')
            line = row[2].strip()+'\t'+row[1].strip()+'\t'+filename[0].strip()+'\t'+row[4].strip()+'\t'+str(row[5]).strip()+'\t'+str(url[0]).strip()+'\t'+md5[0].strip()+'\n'+row[2].strip()+'\t'+row[1].strip()+'\t'+filename[1].strip()+'\t'+row[4].strip()+'\t'+str(row[5]).strip()+'\t'+str(url[1]).strip()+'\t'+md5[1].strip()+'\n'
            with open(output_file, 'a') as writer:
                writer.write(line)
            writer.close()
    return results    
            
def fix_short_read_record(db,host,port,user,output_file,taxon_id,task,task_list):
    #sql = "select * from short_read_data where species_id = 559292"
    sql = "select ID,url,md5sum,species_id from short_read_data where url not like '%;%'"
    info = fetch_data(sql,db,host,int(port),user,'')
    con = pymysql.connect(
            host=host, user='ensadmin', passwd='ensembl', port=int(port), db=db.strip(), charset='utf8mb4', cursorclass=pymysql.cursors.DictCursor
    )
    results = 0
    url = ''
    md5 = ''
    for row in info:
        results += 1
        sq = "select url,md5sum from rnaseq_data_extra where ID = '"+row[0]+"'"
        info1 = fetch_data(sq,db,host,int(port),user,'')
        for data in info1:
            url = data[0]
            md5 = data[1]
       # new_filename = row[3]+';'+row[3].replace('_1','_2') 
        try:
            with con.cursor() as cur:
                cur.execute('UPDATE short_read_data SET url = %s WHERE species_id = %s AND url = %s' ,  (url, row[3], row[1]))
                cur.execute('UPDATE short_read_data SET md5sum = %s WHERE species_id = %s AND md5sum = %s' ,  (md5, row[3], row[2]))
               # if re.search(r"_2",row[3]):
                #    cur.execute('DELETE from short_read_data WHERE filename = %s AND species_id = %s',  (row[3],row[0]))
                #else:
                 #   cur.execute('UPDATE short_read_data SET filename = %s WHERE species_id = %s AND filename = %s' ,  (new_filename, row[0], row[3]))
        except Exception as error:
            print ('Oops! Genebuild update did not complete successfully. Kindly check the error below\n'+str(error))
            break
        else:
            con.commit()
    print('Numer of records fixed = '+str(results))

def get_candidate_assembly(db,host,port,user,output_file):
    #query to fetch records of all existing assemblies
    sql = "SELECT taxonomy, assembly.assembly_id, chain, version, clade, contig_N50, assembly_level, assembly_name, REPLACE(TRIM(LCASE(species_name)), ' ', '_') AS species FROM assembly JOIN meta USING (assembly_id) WHERE genome_rep = 'full'  AND contig_N50 > 100000 AND (rnaseq_data in ('available','red','amber','green')) AND round((30*total_length)/100) > total_gap_length AND is_current = 1"

    info = fetch_data(sql,db,host,int(port),user,'')
    unique_assemblies_per_species = {} #Dictionary to hold one assembly per species where multiple assemblies exists for a species
    for row in info:
        # Replace whitespace with underscore in assembly name
        assembly_name = re.sub(r"\s+", '_', row[7])

        #Get genus id for species being processed
        genus_id = get_taxonomy_group(row[0])
        
        #Write all initial candidate assemblies retrieved from the database
        with open(output_file, 'a') as writer:
           writer.write(str(row[0])+'\t'+str(row[1])+'\t'+row[2]+'\t'+str(row[3])+'\t'+row[4]+'\t'+str(row[5])+'\t'+row[6]+'\t'+assembly_name+'\t'+row[8]+'\t'+genus_id+'\n')
        writer.close()

        #Store one assembly record per species. This would be used to avoid multiple fastq downloads per species
        unique_assemblies_per_species[row[0]] = str(row[0])+'\t'+str(row[1])+'\t'+row[2]+'\t'+str(row[3])+'\t'+row[4]+'\t'+str(row[5])+'\t'+row[6]+'\t'+assembly_name+'\t'+row[8]+'\t'+genus_id

    #Write all unique assemblies per species to file
    unique_assemblies = output_file + '_uniq'
    for assembly in unique_assemblies_per_species:
        with open(unique_assemblies, 'a') as writer:
           writer.write(unique_assemblies_per_species[assembly]+'\n')
        writer.close()        

def get_taxonomy_group(taxon_id):
    """ Function get to genus id of species  """
    while True:
        try:
            handle2 = Entrez.efetch(db='taxonomy', id=taxon_id, retmode='xml')
            break
        except urllib.error.HTTPError as e:
            if e.code == 400:
                to = int(e.hdrs.get("retry-after", 30))
                print ("Got 400. Retrying after {0:d} seconds.".format(to))
                time.sleep(to)
                continue
            else: raise
    record2 = Entrez.read(handle2, validate=False)
    handle2.close()
    report = ''
    #If the rank is not species, we  continue until the genus id is obtained. This is useful when fetching transcriptomic data at a higher level 
    if record2[0]['Rank'] != 'species':
        #This assumes the parent is now at species level
        get_taxonomy_group(record2[0]['ParentTaxId'])
    else:
        report = record2[0]['ParentTaxId']
    return report

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


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--port', help='Port number for host', required=True)
    parser.add_argument('--user', help='Mysql user', required=True)
    parser.add_argument('--dbname', help='Database to be backed up', required=True)
    parser.add_argument('--server', help='Host server for database', required=True)
    parser.add_argument('--asm_file', help='File name to hold short read data', required=True)
    parser.add_argument('--taxon_id', help='Taxon_id to fetch data for', required=True)
    parser.add_argument('--task_dsc', help='Name of task to perform', required=True)
    args = parser.parse_args()
    host = args.server
    db = args.dbname
    user = args.user
    port = args.port
    output_file = args.asm_file
    taxon_id = args.taxon_id
    task = args.task_dsc
    approved_tasks = {
      'get_short_read_at_species_level': 'Returns a tab-delimited of all annotation ready short read data at species level for the taxon_id given',
      'get_short_read_at_genus_level': 'Returns a tab-delimited of all annotation ready short read data at genus level for the taxon_id given',
      'get_all_short_reads': 'Returns a tab-delimited of all annotation ready short read data for the taxon_id given',
      'get_long_read_at_species_level': 'Returns a tab-delimited of all annotation ready long read data at species level for the taxon_id given',
      'get_long_read_at_genus_level': 'Returns a tab-delimited of all annotation ready long read data at genus level for the taxon_id given',
      'get_all_long_reads': 'Returns a tab-delimited of all annotation ready long read data for the taxon_id given'
    }
    #fix_short_read_record(db,host,port,user,output_file,taxon_id,task,approved_tasks)
    records = generate_report(db,host,port,user,output_file,taxon_id,task,approved_tasks)
    print(str(records) + ' entries exists as ready for use in annotation')
