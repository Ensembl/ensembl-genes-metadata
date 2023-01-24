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

""" Module to run the custom assembly registry process.  """

import os
import errno
import subprocess
from datetime import datetime
from sys import exit
import getpass
import re
import time
from sys import setrecursionlimit
import shutil
import random
import string
import requests
import sqlalchemy as db
from pathlib import Path
import argparse
from Bio import Entrez
Entrez.email = 'ensembl-genebuild@ebi.ac.uk'
import sys
import logging
import urllib.request,json
from collections import OrderedDict
from urllib.error import HTTPError
# sqlalchemy requires MySQLdb but MySQLdb doesn't support Python 3.x
# pymysql can be imported and used instead
import pymysql
pymysql.install_as_MySQLdb()
import eHive

class Register_Asm(eHive.BaseRunnable):
    """ Runnable that registers new assemblies in the assembly registry database. """

    def param_defaults(self):
        """ It sets the parameters default values. """
        return {
            'ena_base_url' : 'https://www.ebi.ac.uk/ena/portal/api/search?display=report',
            'files_domain' : 'domain=read&result=read_run',
            'files_fields' : 'run_accession,study_accession,experiment_accession,sample_accession,secondary_sample_accession,instrument_platform,\
            instrument_model,library_layout,library_strategy,nominal_length,read_count,base_count,fastq_ftp,fastq_aspera,submitted_ftp,fastq_md5,\
            library_source,library_selection,center_name,study_alias,experiment_alias,experiment_title,study_title',
            'sample_domain' : 'domain=sample&result=sample',
            'sample_fields' : 'accession,secondary_sample_accession,bio_material,cell_line,cell_type,collected_by,collection_date,country,cultivar,\
            culture_collection,description,dev_stage,ecotype,environmental_sample,first_public,germline,identified_by,isolate,isolation_source,\
            location,mating_type,serotype,serovar,sex,submitted_sex,specimen_voucher,strain,sub_species,sub_strain,tissue_lib,tissue_type,variety,\
            tax_id,scientific_name,sample_alias,center_name,protocol_label,project_name,investigation_type,experimental_factor,sample_collection,sequencing_method',
            'download_method' : 'ftp',
            'separator' : '\t',
            '_read_length' : 1,
            '_centre_name' : 'ENA',
            'print_all_info' : 0,
            'paired_end_only' : 1,
            'ftp_base_dir ' : '/genomes/all/',
            'ftphost' : 'ftp.ncbi.nlm.nih.gov',
            'ftpuser' : 'anonymous',
            'ftppassword' : ''
        }

    def fetch_input(self):
        """ It fetches the input parameters and checks that they are correct. """
        cfg_path = Path(self.param('config_file'))
        registry_settings = {} 
        metazoa_import_types = ['import_community','import_flybase','import_genbank','import_refseq','import_veupathdb','import_wormbase']
        #make sure config file exists before attempting to process it
        if cfg_path.is_file():
            with open(cfg_path) as file:
                while True:
                    lines = file.readline()
                    if not lines:
                        break
                    lines = lines.split('=')
                    registry_settings[lines[0].strip()] = (lines[1].replace('\n','')).strip()
            if metazoa_import_types.count(str(registry_settings['import_type'])) > 0:
                print('Metazoa registration selected')
            else:
                print('Genebuild registration selected')
                registry_settings['import_type'] = 'unannotated'
            self.param('registry', registry_settings)
            self.param('metazoa_import_types', metazoa_import_types)

    #Method to write to database
    def store_assembly(self, database, host, port, user, password, records):
        con = pymysql.connect(
                host=host, user=user, passwd=password, port=port, database=database.strip(), charset='utf8mb4', cursorclass=pymysql.cursors.DictCursor
            )
        chain_version = records[1].split('.')
        new_assemblies = self.param('na')
        updated_assemblies = self.param('ua')
        assemblies_existing_species = self.param('ea')
        n = 0
        u = 0
        e = 0
        if int(records[27][1]) == 1:
            #Store assembly record for the first time into the meta database
            try:
                with con.cursor() as cur:
                    cur.execute('INSERT INTO species_space_log  VALUES (%s,%s,%s)',  (records[11], records[27][0], records[7]))
                    cur.execute('INSERT INTO assembly (chain, version, stable_id_space_id, species_prefix, pri_asm_group, clade, species_id, taxonomy, assembly_path, \
                                genome_rep, rnaseq_data) values (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)',  (chain_version[0], chain_version[1], records[27][0], \
                                records[26], records[24], records[30], records[31], records[11], records[8], records[10], records[23]))
            except Exception as error:
                raise eHive.CompleteEarlyException('Looks like something did not go down well. Kindly check below\n'+str(error))
            else:
                con.commit()
            sql = "SELECT assembly_id FROM assembly WHERE CONCAT(chain,'.',version) = '" + records[1] + "'"
            assembly_id = self.fetch_data(sql,database,host,port,user,password)
            if assembly_id:
                try:
                    with con.cursor() as cur:
                        cur.execute('INSERT INTO meta (assembly_id, assembly_name, assembly_level, wgs_id, scaffold_N50, contig_N50, contig_count, scaffold_count, \
                                    assembly_type, refseq_accession, assembly_date, top_level_count, total_length, total_gap_length, species_name,common_name,\
                                    sec_asm_group, submitter) values (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)', (assembly_id, records[4], records[6],\
                                    records[2], records[16], records[12], records[13], records[15], records[5], records[17], records[18], records[21],  records[19], \
                                    records[20],  records[7],  records[22],  records[25],  records[9]))
                except Exception as error:
                    #Rollback changes already made to assembly and species space tables
                    if int(records[27][1]) > 1:
                        space_id = records[27][0] - 1
                        cur.execute('UPDATE species_space_log SET current_space = %s WHERE species_id = %s',  (space_id, records[11]))
                        cur.execute('DELETE from assembly WHERE assembly_id = %s',  (assembly_id))
                        con.commit()
                    else:
                        cur.execute('DELETE from species_space_log WHERE species_id = %s',  (records[11]))
                        cur.execute('DELETE from assembly WHERE assembly_id = %s',  (assembly_id))
                        con.commit()
                    raise eHive.CompleteEarlyException('Looks like something did not go down well. Kindly check below\n'+str(error))
                else:
                    con.commit()
            #If we get here, it means new assembly registration completed without error
            #Create a list of all newly registered assemblies for reporting purposes via Slack
            new_assemblies[records[1]] = records[1] + '\t' + records[30] + '\t' + records[7] + '\t' + records[22] + '\t' + records[6] + '\t' + records[4] + '\t' + records[9]
            self.param('na',new_assemblies)
        elif int(records[27][1]) == 2:
            try:
                with con.cursor() as cur:
                    cur.execute('INSERT INTO assembly (chain, version, stable_id_space_id, species_prefix, pri_asm_group, clade, species_id, taxonomy, assembly_path,\
                                genome_rep, rnaseq_data) values (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)',  (chain_version[0], chain_version[1], records[27][0], records[26],\
                                records[24],records[30],records[31], records[11], records[8], records[10], records[23]))
            except Exception as error:
                raise eHive.CompleteEarlyException('Looks like something did not go down well. Kindly check below\n'+str(error))
            else:
                con.commit()
            sql = "SELECT assembly_id FROM assembly WHERE CONCAT(chain,'.',version) = '" + records[1] + "'"
            assembly_id = self.fetch_data(sql,database,host,port,user,password)
            if assembly_id:
                try:
                    with con.cursor() as cur:
                        cur.execute('INSERT INTO meta (assembly_id, assembly_name, assembly_level, wgs_id, scaffold_N50, contig_N50, contig_count, scaffold_count, \
                                    assembly_type, refseq_accession, assembly_date, top_level_count, total_length, total_gap_length, species_name,common_name, \
                                    sec_asm_group, submitter) values (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)', (assembly_id, records[4], records[6], \
                                    records[2], records[16], records[12], records[13], records[15], records[5], records[17], records[18], records[21],  records[19], \
                                    records[20],  records[7],  records[22],  records[25],  records[9]))
                except Exception as error:
                    #Rollback changes already made to assembly and species space tables
                    cur.execute('DELETE from assembly WHERE assembly_id = %s',  (assembly_id))
                    con.commit()
                    raise eHive.CompleteEarlyException('Looks like something did not go down well. Kindly check below\n'+str(error))
                else:
                    con.commit()
            #If we get here, it means updated assembly registration completed without error
            #Create a list of all updated assemblies for reporting purposes via Slack
            updated_assemblies[records[1]] = records[1] + '\t' + records[30] + '\t' + records[7] + '\t' + records[22] + '\t' + records[6] + '\t' + records[4] + '\t' + chain_version[1] + '\t' + records[9]
            self.param('ua',updated_assemblies)
        #Updating existing assembly with a newer version whilst maintaining same stable space range
        elif int(records[27][1]) == 3:
            try:
                with con.cursor() as cur:
                    cur.execute('INSERT INTO stable_id_space  VALUES (%s,%s,%s)',  (records[27][0], records[28], records[29]))
                    cur.execute('UPDATE species_space_log SET current_space = %s WHERE species_id = %s',  (records[27][0], records[11]))
                    cur.execute('INSERT INTO assembly (chain, version, stable_id_space_id, species_prefix, pri_asm_group, clade, species_id, taxonomy, assembly_path, \
                                genome_rep,rnaseq_data) values (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)', (chain_version[0], chain_version[1], records[27][0], records[26],\
                                records[24], records[30],records[31], records[11], records[8], records[10], records[23]))
            except Exception as error:
                raise eHive.CompleteEarlyException('Looks like something did not go down well. Kindly check below\n'+str(error))
            else:
                con.commit()
            #Use assembly id for newly registered assembly to record assembly's metadata
            sql = "SELECT assembly_id FROM assembly WHERE CONCAT(chain,'.',version) = '" + records[1] + "'"
            assembly_id = self.fetch_data(sql,database,host,port,user,password)
            if assembly_id:
                try:
                    with con.cursor() as cur:
                        cur.execute('INSERT INTO meta (assembly_id, assembly_name, assembly_level, wgs_id, scaffold_N50, contig_N50, contig_count, scaffold_count,\
                                    assembly_type, refseq_accession, assembly_date, top_level_count, total_length, total_gap_length, species_name,common_name, \
                                    sec_asm_group, submitter) values (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)', (assembly_id, records[4], records[6],\
                                    records[2], records[16], records[12], records[13], records[15], records[5], records[17], records[18], records[21],  records[19], \
                                    records[20],  records[7],  records[22],  records[25],  records[9]))
                except Exception as error:
                     #Rollback changes already made to assembly and species space tables
                     space_id = records[27][0] - 1
                     cur.execute('UPDATE species_space_log SET current_space = %s WHERE species_id = %s',  (space_id, records[11]))
                     cur.execute('DELETE from assembly WHERE assembly_id = %s',  (assembly_id))
                     con.commit()
                     raise eHive.CompleteEarlyException('Looks like something did not go down well. Kindly check below\n'+str(error))
                else:
                    con.commit()
            #If we get here, it means registration for new assemblies of existing species completed without error
            #Create a list of all registered assemblies for reporting purposes via Slack
            assemblies_existing_species[records[1]] = records[1] + '\t' + records[30] + '\t' + records[7] + '\t' + records[22] + '\t' + records[6] + '\t' + records[4] + '\t' + records[9]
            self.param('ea',assemblies_existing_species)
        #Registering a new assembly on a different chain for an existing species. A new stable space range is required here
        elif int(records[27][1]) == 4:
            #Re-use existing stable space range where assembly is of different chain from existing ones belonging to same species
            try:
                with con.cursor() as cur:
                    cur.execute('UPDATE species_space_log SET current_space = %s WHERE species_id = %s',  (records[27][0], records[11]))
                    cur.execute('INSERT INTO assembly (chain, version, stable_id_space_id, species_prefix, pri_asm_group, clade, species_id, taxonomy, assembly_path, \
                                genome_rep, rnaseq_data) values (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)',  (chain_version[0], chain_version[1], records[27][0], records[26], \
                                records[24], records[30],records[31], records[11], records[8], records[10], records[23]))
            except Exception as error:
                raise eHive.CompleteEarlyException('Looks like something did not go down well. Kindly check below\n'+str(error))
            else:
                con.commit()
            #Use assembly id for newly registered assembly to record assembly's metadata
            sql = "SELECT assembly_id FROM assembly WHERE CONCAT(chain,'.',version) = '" + records[1] + "'"
            assembly_id = self.fetch_data(sql,database,host,port,user,password)
            if assembly_id:
                try:
                    with con.cursor() as cur:
                        cur.execute('INSERT INTO meta (assembly_id, assembly_name, assembly_level, wgs_id, scaffold_N50, contig_N50, contig_count, scaffold_count, \
                                    assembly_type, refseq_accession, assembly_date, top_level_count, total_length, total_gap_length, species_name,common_name, \
                                    sec_asm_group, submitter) values (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)', (assembly_id, records[4], records[6], \
                                    records[2], records[16], records[12], records[13], records[15], records[5], records[17], records[18], records[21],  records[19], \
                                    records[20],  records[7],  records[22],  records[25],  records[9]))
                except Exception as error:
                    #Rollback changes already made to assembly and species space tables
                    space_id = records[27][0] - 1
                    cur.execute('UPDATE species_space_log SET current_space = %s WHERE species_id = %s',  (space_id, records[11]))
                    cur.execute('DELETE from assembly WHERE assembly_id = %s',  (assembly_id))
                    con.commit()
                    raise eHive.CompleteEarlyException('Looks like something did not go down well. Kindly check below\n'+str(error))
                else:
                    con.commit()
            #If we get here, it means registration for new assemblies of existing species completed without error
            #Create a list of all registered assemblies for reporting purposes via Slack
            assemblies_existing_species[records[1]] = records[1] + '\t' + records[30] + '\t' + records[7] + '\t' + records[22] + '\t' + records[6] + '\t' + records[4] + '\t' + records[9]
            self.param('ea',assemblies_existing_species)
        return [new_assemblies,updated_assemblies,assemblies_existing_species,'0']
    def slack_reporting(self,new,updates,diff,imp):
        #Reporting assembly registration in categories
        if len(imp) > 0:
            cnt = 0
            asm_meta = ''
            report = "Accession\tImporter\tAnnotation method\n"
            for entry in new.values():
                cnt+=1
                asm_meta = asm_meta + entry + '\n'
            report = 'Total number of imported annotations = ' + str(cnt) + '\n' + report + str(asm_meta)
            payload="{\"channel\": \"#genebuildregistry\", \"username\": \"registry_messenger\", \"text\": \"" + report  +"\"}"
            url = os.getenv('slack_token')
            headers = {'content-type': 'application/json', 'Accept-Charset': 'UTF-8'}
            r = requests.post(url, data=payload, headers=headers)
        elif len(new) > 0:
            cnt = 0
            asm_meta = ''
            report = "Accession\tClade\tSpecies name\tCommon name\tAssembly level\tAssembly name\tSubmitter\n"
            for entry in new.values():
                cnt+=1
                asm_meta = asm_meta + entry + '\n'
            report = 'Total number of new species registered = ' + str(cnt) + '\n' + report + str(asm_meta)
            payload="{\"channel\": \"@denye\", \"username\": \"registry_messenger\", \"text\": \"" + report  +"\"}"
            url = os.getenv('slack_token')
            headers = {'content-type': 'application/json', 'Accept-Charset': 'UTF-8'}
            r = requests.post(url, data=payload, headers=headers)
        if len(updates) > 0:
            cnt = 0
            asm_meta = ''
            report = "Accession\tClade\tSpecies name\tCommon name\tAssembly level\tAssembly name\tVersion\tSubmitter\n"
            for entry in updates.values():
                cnt+=1
                asm_meta = asm_meta + entry + '\n'
            report = 'Total number of assembly updates on existing species  = ' + str(cnt) + '\n' + report + str(asm_meta)
            payload="{\"channel\": \"@denye\", \"username\": \"registry_messenger\", \"text\": \"" + report  +"\"}"
            url = os.getenv('slack_token')
            headers = {'content-type': 'application/json', 'Accept-Charset': 'UTF-8'}
            r = requests.post(url, data=payload, headers=headers)
        if len(diff) > 0:
            cnt = 0
            asm_meta = ''
            report = "Accession\tClade\tSpecies name\tCommon name\tAssembly level\tAssembly name\tSubmitter\n"
            for entry in diff.values():
                cnt+=1
                asm_meta = asm_meta + entry + '\n'
            report = 'Total number of new assemblies on existing species = ' + str(cnt) + '\n' + report + str(asm_meta)
            payload="{\"channel\": \"@denye\", \"username\": \"registry_messenger\", \"text\": \"" + report  +"\"}"
            url = os.getenv('slack_token')
            headers = {'content-type': 'application/json', 'Accept-Charset': 'UTF-8'}
            r = requests.post(url, data=payload, headers=headers)

    #Method to establish connection to database and query it
    def fetch_data(self, query, database, host, port, user, password):
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

    def run(self):
        """ Perform actual registration. """
        registry_settings = self.param('registry')
        #query to fetch records of all existing assemblies
        sql = "SELECT species_prefix, concat(chain,'.', version), taxonomy FROM assembly order by chain,version"
        info = self.fetch_data(sql,registry_settings['registry_dbname'],registry_settings['registry_db_host'],int(registry_settings['registry_db_port']),\
                               registry_settings['user_w'],registry_settings['password'])
        existing_gca = {}
        existing_prefix = {}
        existing_taxons = {}
        existing_annotations = {}
        new_assemblies = {}
        updated_assemblies = {}
        existing_assemblies = {}
        self.param('na',new_assemblies)
        self.param('ua',updated_assemblies)
        self.param('ea',existing_assemblies)
        pattern = ['[',']']
        res = []
        metazoa_import_types = self.param('metazoa_import_types')
        #convert input string of accessions into a list
        acc_list = registry_settings['assembly_accessions'].translate( {ord(elem): None for elem in pattern} ).split(',')
        acc_format = re.compile(r"^GC[AF]_\d+\.\d+$") #regex to filter genbank accession
        string = ''.join(acc_list) #convert list of accessions to a string to enable pattern search
        match_list = re.findall(acc_format, string)
        #Check all accessions given are in right format
        for i in range(0, len(acc_list)):
            if not re.match(acc_format, acc_list[i]):
                raise ValueError('Accession '+ acc_list[i] + ' not in the right format ')
        #Create a dictionary of existing assemblies and species prefixes. This avoids duplicates in the database and registration process
        for row in info:
            existing_gca[row[1]] = row[0]
            existing_taxons[str(row[2])] = row[0]
            existing_prefix[row[0]] = row[2]
        ecnt = len(existing_prefix)
        self.param('existing_prefix',existing_prefix)
        sql = "SELECT assembly_accession, annotation_source FROM genebuild_status"
        annotations = self.fetch_data(sql,registry_settings['registry_dbname'],registry_settings['registry_db_host'],int(registry_settings['registry_db_port']),\
                                      registry_settings['user_w'],registry_settings['password'])
        #Create a disctionary of existing annotations.
        for row in annotations:
            existing_annotations[row[0]] = row[1]

        #Retrieve current stable space assignment per species
        sql = "SELECT species_id, current_space FROM species_space_log"
        info = self.fetch_data(sql,registry_settings['registry_dbname'],registry_settings['registry_db_host'],int(registry_settings['registry_db_port']),\
                               registry_settings['user_w'],registry_settings['password'])
        existing_stable_id = {}
        for space_id in info:
            existing_stable_id[str(space_id[0])] = space_id[1]

        #Retrieve current stable space id range
        sql = "SELECT stable_id_space_id, stable_id_space_start, stable_id_space_end FROM stable_id_space"
        info = self.fetch_data(sql,registry_settings['registry_dbname'],registry_settings['registry_db_host'],int(registry_settings['registry_db_port']),\
                               registry_settings['user_w'],registry_settings['password'])
        existing_stable_space_range = {}
        for space_range in info:
            existing_stable_space_range[space_range[0]] = str(space_range[1]) + '\t' + str(space_range[2])
        
        self.param('stable_id_space_range',existing_stable_space_range)
        self.param('species_space_id',existing_stable_id)
        self.param('existing_species',existing_taxons)
        registry_update = []
        assembly_id = {}
        #Check if metazoa performing custom assembly registration or genebuild
        if registry_settings['import_type'] in metazoa_import_types:
            #Get all existing assembly ids
            sql = "SELECT CONCAT(chain,'.',version), assembly_id FROM assembly"
            results = self.fetch_data(sql,registry_settings['registry_dbname'],registry_settings['registry_db_host'],int(registry_settings['registry_db_port']),\
                                      registry_settings['user_w'],registry_settings['password'])
            for row in results:
                assembly_id[row[0]] = row[1]

            #Get a list of new assemblies to register and check against existing entries
            new_record = 0
            imported_assemblies = {}
            for accession in acc_list:
                if (accession in existing_gca):
                    #Check that assembly with same annotation source does exist in the annotation/genebuild status table
                    sql = "SELECT DISTINCT(annotation_source) FROM genebuild_status WHERE assembly_accession = '" + accession + "'"
                    results = self.fetch_data(sql,registry_settings['registry_dbname'],registry_settings['registry_db_host'],int(registry_settings['registry_db_port']),\
                                              registry_settings['user_w'],registry_settings['password'])
                    if registry_settings['import_type'] in str(results):
                        print ('Annotation for assembly '+ accession + ' already exists. No further action is needed.')
                        continue
                    else:
                        #Retrieve all annotation records for the assembly from the genebuild_status table
                        sql = "SELECT * FROM genebuild_status WHERE assembly_accession = '" + accession + "'"
                        annotation_record = self.fetch_data(sql,registry_settings['registry_dbname'],registry_settings['registry_db_host'],\
                                                            int(registry_settings['registry_db_port']),registry_settings['user_w'],registry_settings['password'])
                        if annotation_record:
                            print('The following annotation(s) exists for this assembly')
                            for record in annotation_record:
                                annotation_record = record
                                print('Accession: ' + accession + '\tStatus: ' + annotation_record[1] + '\tDate started: '+str(annotation_record[2]) +\
                                      '\tDate completed: '+str(annotation_record[3]) + '\tGenebuiler: '+str(annotation_record[4]) + '\tAnnotation method: '+\
                                      str(annotation_record[8]))
                            try:
                                response = input('Do you wish to continue with the registration of this assembly '+accession+'? Y/N\n')
                            except Exception as EOFError:
                                raise eHive.CompleteEarlyException('No response was received. Registration terminated.\n')
                            if (response.isalpha() and response.upper() == 'Y'):
                                #Add assembly import to the annotation status table
                                con = pymysql.connect(
                                                      host=registry_settings['registry_db_host'], user=registry_settings['user_w'], passwd=registry_settings['password'], \
                                                      port=int(registry_settings['registry_db_port']), database=registry_settings['registry_dbname'].strip(), charset='utf8mb4',\
                                                      cursorclass=pymysql.cursors.DictCursor
                                )
                                try:
                                    with con.cursor() as cur:
                                        cur.execute('INSERT INTO genebuild_status (assembly_accession, progress_status, date_started, date_completed, genebuilder, \
                                                     assembly_id,is_current, annotation_source) VALUES (%s,%s,%s,%s,%s,%s,%s,%s)',(accession, 'completed', \
                                                     datetime.today().strftime('%Y-%m-%d'),datetime.today().strftime('%Y-%m-%d'), getpass.getuser(), assembly_id[accession],\
                                                     1, registry_settings['import_type']))
                                except Exception as error:
                                    raise eHive.CompleteEarlyException('Looks like something did not go down well. Kindly check below\n'+str(error))
                                else:
                                    con.commit()
                                    print('Genome annotation for assembly '+accession+' successfully imported into registry')
                                    #Create a list of all newly registered assemblies for reporting purposes via Slack
                                    imported_assemblies[accession] = accession + '\t' + getpass.getuser() + '\t' + registry_settings['import_type']
                            else:
                                 continue
                        else:
                            #Add assembly import to the annotation status table for the first time as no entry exists for the assembly yet
                            print('No genebuild entry for assembly '+str(accession))
                            #Add assembly import to the annotation status table
                            con = pymysql.connect(
                                                  host=registry_settings['registry_db_host'], user=registry_settings['user_w'], passwd=registry_settings['password'], \
                                                  port=int(registry_settings['registry_db_port']), database=registry_settings['registry_dbname'].strip(), charset='utf8mb4',\
                                                  cursorclass=pymysql.cursors.DictCursor)
                            try:
                                with con.cursor() as cur:
                                    cur.execute('INSERT INTO genebuild_status (assembly_accession, progress_status, date_started, date_completed, genebuilder, assembly_id,\
                                                is_current, annotation_source) VALUES (%s,%s,%s,%s,%s,%s,%s,%s)',(accession, 'completed', datetime.today().strftime('%Y-%m-%d'),\
                                                datetime.today().strftime('%Y-%m-%d'), getpass.getuser(), assembly_id[accession], 1, registry_settings['import_type']))
                            except Exception as error:
                                raise eHive.CompleteEarlyException('Error updating the genebuild status of this assembly. Kindly check below\n'+str(error))
                            else:
                                con.commit()
                                print('Genome annotation for assembly '+accession+' successfully imported into registry')
                                imported_assemblies[accession] = accession + '\t' + getpass.getuser() + '\t' + registry_settings['import_type']
                        res = [imported_assemblies,updated_assemblies,existing_assemblies,'1']
                else:
                    #Register assembly as standard since it does not already exist. Also add it to the annotation status table
                    print ('Registration in progress for assembly ' + accession)
                    self.register_assembly(accession)
                    res = self.store_assembly(registry_settings['registry_dbname'],registry_settings['registry_db_host'],int(registry_settings['registry_db_port']),\
                                                   registry_settings['user_w'],registry_settings['password'],self.param('metadata'))
                    #After registering new assembly for metazoa, it is necessary to update the list of existing assemblies before attempting to register it's genebuild status
                    #Get all existing assembly ids
                    sql = "SELECT assembly_id FROM assembly WHERE CONCAT(chain,'.',version) = '"+accession+"'"
                    results = self.fetch_data(sql,registry_settings['registry_dbname'],registry_settings['registry_db_host'],int(registry_settings['registry_db_port']),\
                                      registry_settings['user_w'],registry_settings['password'])
                    for row in results:
                        assembly_id[accession] = row[0]
                    #Add assembly import to the annotation status table
                    con = pymysql.connect(
                                          host=registry_settings['registry_db_host'], user=registry_settings['user_w'], passwd=registry_settings['password'], \
                                          port=int(registry_settings['registry_db_port']), database=registry_settings['registry_dbname'].strip(), charset='utf8mb4',\
                                          cursorclass=pymysql.cursors.DictCursor)
                    try:
                        with con.cursor() as cur:
                            cur.execute('INSERT INTO genebuild_status (assembly_accession, progress_status, date_started, date_completed, genebuilder, assembly_id,\
                                        is_current, annotation_source) VALUES (%s,%s,%s,%s,%s,%s,%s,%s)',(accession, 'completed', datetime.today().strftime('%Y-%m-%d'),\
                                        datetime.today().strftime('%Y-%m-%d'), getpass.getuser(), assembly_id[accession], 1, registry_settings['import_type']))
                    except Exception as error:
                        raise eHive.CompleteEarlyException('Error updating the genebuild status of this assembly. Kindly check below\n'+str(error))
                    else:
                        con.commit()
                        print('Genome annotation for assembly '+accession+' successfully imported into registry')
                        imported_assemblies[accession] = accession + '\t' + getpass.getuser() + '\t' + registry_settings['import_type']
                        res = [imported_assemblies,updated_assemblies,existing_assemblies,'1']
                    registry_update.append(self.param('metadata'))
                    existing_gca[accession] = self.param('metadata')[26]
                    existing_taxons[self.param('metadata')[11]] = self.param('metadata')[26]
                    self.param('existing_species',existing_taxons)
                  
        else: 
            #Registration being performed by genebuild
            res = []
            #Get a list of new assemblies to register and check against existing entries
            for accession in acc_list:
                if (accession in existing_gca):
                    print ('Assembly '+ accession + ' already exists. No further action is needed.')
                    continue
                else:
                    print ('Registration in progress for assembly ' + accession)
                    self.register_assembly(accession)
                    #Temporarily write to db per record instead of originally designed to write in bulk
                    res = self.store_assembly(registry_settings['registry_dbname'],registry_settings['registry_db_host'],int(registry_settings['registry_db_port']),\
                                registry_settings['user_w'],registry_settings['password'],self.param('metadata'))
                    registry_update.append(self.param('metadata'))
                    print('Updating existing accessions to include '+str(accession)+' and '+str(self.param('metadata')[26]))
                    existing_gca[accession] = self.param('metadata')[26]
                    existing_taxons[self.param('metadata')[11]] = self.param('metadata')[26]
                    self.param('existing_species',existing_taxons)
        self.param('records', registry_update)
        if res:
            self.slack_reporting(res[0],res[1],res[2],res[3])
        else:
            raise eHive.CompleteEarlyException('No new assembly or update at this time.\n')
    
    #For linked assemblies, find the id associated with the assembly of interest
    def get_ids(self,term):
        ids = []
        handle = Entrez.esearch(db="assembly", term=term)
        record = Entrez.read(handle)
        ids.append(record["IdList"])
        return ids

    #Method to process assembly registration 
    def register_assembly(self,accession):
        """ `function that takes genbank accessions and finds all required metadata to store in registry database """
        #Query public archives for assembly related information as per gca supplied
        registry_settings = self.param('registry')
        existing_taxons = self.param('existing_species')
        existing_space_range = self.param('stable_id_space_range')
        asm_group_dict = {'genomes25' : 'ebp', 'ebp' : 'None', 'abp' : 'ebp', 'asg' : 'ebp', 'b10k' : 'ebp', 'dtol' : 'erga ebp', 'erga' : 'ebp', 'ergapp' : 'erga', \
                          'g10k' : 'erga  ebp', 'tol': 'None', 'ungrouped' : 'None', 'vgp' : 'ebp'} 
        for id in self.get_ids(accession):
            if len(id) > 1: #Assembly is linked to other assembly(ies)
                for i in range(len(id)): #Get metadata for each linked assembly
                    assembly_metadata = self.get_assembly_metadata(id[i],accession)
                    if assembly_metadata[0]: #stop searching linked assemblies once matching metadata is found
                        break
            else:
                #Assembly is not linked to any other assembly
                assembly_metadata = self.get_assembly_metadata(id,accession)
        #Get clade related information for species
        clade_rep = self.get_taxonomy_group(assembly_metadata[11],assembly_metadata[8],assembly_metadata[4],accession)
        species_details =  clade_rep.split(':') #Retrieve taxonomy related data 
        if species_details[1]:
            common_name = species_details[1]
        else:
            #This means all attempts have been made to find a meaningful common name for this species and in the absence of none, we assign NA
            common_name = 'NA'
        #Check transcriptomic data status at species level
        trans_status = self.transcriptomic_status(species_details[3])
        #Check what major project assembly falls into for proper grouping
        assembly_group = self.set_assembly_group(assembly_metadata[3])
        if asm_group_dict[assembly_group]:
            sec_asm_group = asm_group_dict[assembly_group]
        if existing_taxons.get(str(assembly_metadata[11].strip())) is not None:
            #Use existing prefix for species
            species_prefix= existing_taxons[assembly_metadata[11]]
        else:
            #Generate initial species prefix for new species
            pattern = "[\n|\r|\r\n|\s+|\d+]" #pattern to remove new line character, whitespace and digits from species name
            species_name = re.sub(pattern, '', assembly_metadata[7] )
            # Create a regex pattern to replace all non letters or special characters in string with empty string
            species_name = re.sub("[^A-Z]", "", species_name,0,re.IGNORECASE)
            pattern = r'[' + string.punctuation + ']' #replace punctuation
            species_name = re.sub(pattern, '', species_name )
            #Randomly generate 3 letters to be appended to the default ENS prefix
            species_prefix = ''.join(random.choice(species_name) for i in range(3)).upper()
            species_prefix = 'ENS' + species_prefix   
            #Check if new prefix is unique
            self.generate_species_prefix(species_name,species_prefix,assembly_metadata[11])
            species_prefix = self.param('unique_prefix')
        #Get stable space range for assembly
        stable_space_id = self.get_stable_space_id(assembly_metadata[11],accession).split('\t')
        stable_id_start_end = []
        if str(stable_space_id[1]) == '3':
            #This requires generating a new stable space range
            current_space_id = int(stable_space_id[0]) - 1
            space_range = existing_space_range[current_space_id].split('\t')
            if (assembly_metadata[11] == 7955):
                stable_id_start_end.append(int(space_range[1]) + 1)
                stable_id_start_end.append(int(space_range[1]) + 15000000)
            elif (assembly_metadata[11] == 10181):
                stable_id_start_end.append(int(space_range[1]) + 1)
                stable_id_start_end.append(int(space_range[1]) + 10000000)
            elif (assembly_metadata[11] != 10029):
                stable_id_start_end.append(int(space_range[1]) + 1)
                stable_id_start_end.append(int(space_range[1]) + 10000000)
            else:
                stable_id_start_end.append(int(space_range[1]) + 1)
                stable_id_start_end.append(int(space_range[1]) + 5000000)
            existing_space_range[str(stable_space_id[0])] = str(stable_id_start_end[0]) + '\t' + str(stable_id_start_end[1]) 
        else:
            stable_id_start_end = existing_space_range[int(stable_space_id[0])].split('\t')
        #Add extra fields to the list of meta data values
        existing_stable_id = self.param('species_space_id')
        existing_stable_id[str(assembly_metadata[11])] = int(stable_space_id[0])
        self.param('stable_id_space_range',existing_space_range)
        self.param('species_space_id',existing_stable_id)
        assembly_metadata.extend([common_name,trans_status,assembly_group,sec_asm_group,species_prefix,stable_space_id,stable_id_start_end[0],stable_id_start_end[1],\
                                 species_details[4],species_details[3]])
        self.param('metadata',assembly_metadata)
    
    #Function to generate species prefix
    def generate_species_prefix(self,species_name,species_prefix,taxon_id):
        existing_prefix = self.param('existing_prefix')
        pcnt = len(existing_prefix)
        if species_prefix in existing_prefix.keys():
            #Regenerate another species prefix
            species_prefix = ''.join(random.choice(species_name) for i in range(3)).upper()
            species_prefix = 'ENS' + species_prefix
            self.generate_species_prefix(species_name,species_prefix,taxon_id)
        else:
            #Unique prefix generated
            existing_prefix[str(species_prefix)] = taxon_id
            self.param('existing_prefix',existing_prefix)
            pcnt+=1
            self.param('unique_prefix',species_prefix)
    
    #Function to generate stable space values
    def get_stable_space_id(self,taxon_id,accession):
        registry_settings = self.param('registry')
        existing_stable_space_range = self.param('stable_id_space_range')
        existing_stable_id_assignment = self.param('species_space_id')
        chain = accession.split('.')
        if taxon_id in existing_stable_id_assignment.keys():
            #for entry in existing_stable_id_assignment.keys():
            sql = "SELECT MAX(version), stable_id_space_id FROM assembly WHERE chain = '" + chain[0] + "'"
            info = self.fetch_data(sql,registry_settings['registry_dbname'],registry_settings['registry_db_host'],int(registry_settings['registry_db_port']),\
                                   registry_settings['user_w'],registry_settings['password'])
            #Check if this is a valid update
            for row in info:
                if row[0] is None:
                    print('Registering a different chain '+str(accession))
                    #Assembly has a different chain from other species related assemblies, therefore assign a new stable id space range
                    for key, val in existing_stable_id_assignment.items():
                        if str(key) == str(taxon_id):
                            stable_id_space_id = val + 1 #val + 1
                            if str(stable_id_space_id) in str(existing_stable_space_range.keys()):
                                rgn_flag = 4
                                return str(stable_id_space_id) + '\t' + str(rgn_flag)
                            else:
                                rgn_flag = 3
                                return str(stable_id_space_id) + '\t' + str(rgn_flag)
                elif row[0] >= int(chain[1]):
                    print('Only a new/updated assembly can be registered. Check database for existing assembly information')
                    return 0
                else:
                    print('Re-using same space')
                    rgn_flag = 2
                    #Perform version update by re-using same stable is space range
                    return str(row[1]) + '\t' + str(rgn_flag)
        else:
             #Assembly/species is being registered for the first time
            stable_id_space_id = 1
            rgn_flag = 1
            #Update current space assignment for each new species registered
            existing_stable_id_assignment[str(taxon_id)] = int(stable_id_space_id)
            self.param('species_space_id', existing_stable_id_assignment)
            return str(stable_id_space_id) + '\t' + str(rgn_flag)
 
    #Function to classify assembly into clade
    def get_taxonomy_group(self,taxon_id,assembly_path,assembly_name,accession):
        """ Function to retrieve/classify species into clades based on taxon_id """
        # Now connect NCBI again using the tax_id
        # Entrez.efetch will give you various information
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
        clade = ''
        lineage_rpt = {'Rodentia': 'rodentia','Primates': 'primates','Mammalia': 'mammalia','Amphibia': 'amphibians','Teleostei': 'teleostei','Marsupialia': 'marsupials',\
                       'Aves': 'aves','Sauropsida': 'reptiles','Chondrichthyes': 'sharks','Eukaryota': 'non_vertebrates','Metazoa': 'metazoa','Viral': 'viral',\
                       'Viruses': 'viral','Viridiplantae': 'plants','Arthropoda': 'arthropods','Lepidoptera': 'lepidoptera','Insecta': 'insects','Hymenoptera': 'hymenoptera',\
                       'Hemiptera': 'hemiptera','Coleoptera': 'coleoptera','Diptera': 'diptera','Mollusca': 'mollusca','Vertebrata': 'vertebrates','Alveolata': 'protists',\
                       'Amoebozoa': 'protists','Choanoflagellida':'protists','Cryptophyta':'protists','Euglenozoa':'protists','Fungi':'fungi','Fornicata':'protists',\
                       'Heterolobosea':'protists','Parabasalia':'protists','Rhizaria':'protists','Stramenopiles':'protists','Trichoptera':'trichoptera',\
                       'Plasmodium':'plasmodium','Crustacea':'crustacea','Ascomycota':'ascomycota','Basidiomycota':'basidiomycota'}
        if ('OtherNames' in record2[0]):
            common_name = record2[0]['OtherNames']
            if(('CommonName' in common_name) and (not common_name['CommonName'])):
                report = self.get_common_name_gbif(record2[0]['ScientificName'].replace(' ','%20'))
            elif (('CommonName' in common_name) and (common_name['CommonName'])):
                report = 'common_name:' + common_name['CommonName'][0].capitalize() + ':' + record2[0]['ScientificName']
            elif ('GenbankCommonName' in common_name):
                report = 'common_name:' + common_name['GenbankCommonName'].capitalize() + ':' + record2[0]['ScientificName']
        else:
            report = self.get_common_name_gbif(record2[0]['ScientificName'].replace(' ','%20'))
        #Check if common name was retrieved via gbif, else try using assembly report file
        common_name = report.split(':')
        if common_name[1] == '':
            report = self.get_common_name_via_report_file(assembly_path,assembly_name,accession,record2[0]['ScientificName'])
        #Reverse the dictionary entry returned by API so we have lowere ancestors before higher up ones
        tax_list = record2[0]['LineageEx']
        rev_list = {}
        for key in reversed(tax_list):
            rev_list[key['TaxId']] = key['ScientificName']
        for tax_element in rev_list:
            if rev_list[tax_element] in lineage_rpt.keys():
                clade = lineage_rpt[rev_list[tax_element]]
                break
        
        #If the rank of the species is subspecies, we set the species id to the parent taxon for transcriptomic data selection 
        if record2[0]['Rank'] != 'species': 
            report += ':' + record2[0]['ParentTaxId'] + ':'+clade
        else:
            report += ':' + record2[0]['TaxId'] + ':'+clade
        return report

    
    def get_common_name_via_report_file(self,assembly_path,assembly_name,accession,species_name):
        """ Fetch common name from assembly stats file """
        assembly_name = assembly_name.replace(' ','_')
        url = assembly_path +'/' + accession+'_' + assembly_name + '_assembly_stats.txt'
        report = ""
        for line in urllib.request.urlopen(url):
            if (re.search(r'Organism name',str(line))):
                res = re.findall('\((.*?)\)', line.decode('utf-8'))
                if res[0]:
                    report = 'common_name:' + res[0].capitalize() + ':' + species_name
                    break
                else:
                    report = 'common_name:' + '' + ':' + species_name
        if not report:
            print('Common name could not be determined for this species')
            report = 'common_name:' + '' + ':' + species_name
        return report

    def get_common_name_gbif(self,species_name):
        """ Fetch common name from alternative sources where it cannot be retrieved via NCBI """
        url = "https://api.gbif.org/v1/species/match?name="+species_name
        with urllib.request.urlopen(url) as response:
            data = json.loads(response.read())
            url = "https://api.gbif.org/v1/species/"+str(data['usageKey'])+"/vernacularNames"
        with urllib.request.urlopen(url) as common_name:
            data = json.loads(common_name.read())
            species_name = species_name.replace('%20',' ')
            if ('results' in data):
                res = data['results']
                if res:
                    for entry in res:
                        if (entry['language'] == 'eng' and entry['vernacularName']):
                            report = 'common_name:' + (entry['vernacularName']).capitalize() + ':' + species_name
                            break
                        else:
                            report = 'common_name:' + '' + ':' + species_name
                else:
                    report = 'common_name:' + '' + ':' + species_name
                    print('No common name entry was found for the species from GBIF')
            else:
                 report = 'common_name:' + '' + ':' + species_name
                 print('No common name entry was found for the species from GBIF')
        return report

    def transcriptomic_status(self,taxon_id):
        """ Check status of transcriptomic data both in the registry for exisiting species and publicly for new species """
        total_read = 0
        entry = 0
        status = 'not available'
        ftp_path = self.param('ena_base_url') + '&query="tax_tree(' + taxon_id + ') AND instrument_platform=ILLUMINA AND library_source=TRANSCRIPTOMIC"&' \
                   + self.param('files_domain') + '&fields=' +self.param('files_fields')
        response = requests.get(ftp_path)
        fastq_file = 'fastq_'+self.param('download_method')
        header_pattern = '^(\w+.*)$'
        content_pattern = '^[a-z]'
        encoding_style = response.encoding
        unwanted_reads = []
        wanted_reads = []
        #Check registry for previously set transcriptomic status for species
        registry_settings = self.param('registry')
        #query to fetch records of all existing assemblies
        sql = "SELECT DISTINCT(rnaseq_data) from assembly where species_id = " + taxon_id        
        info = self.fetch_data(sql,registry_settings['registry_dbname'],registry_settings['registry_db_host'],int(registry_settings['registry_db_port']),\
                               registry_settings['user_w'],registry_settings['password'])
        #Check the status of transcriptomic data for this species in the registry if any before attempting to look for new data
        if info:
            for row in info:
                status = row[0]
                if (re.search(r'[red|green|amber|available]',status)):
                    return status
        if response:
            read_file = registry_settings['output_path'] + "/ena_read_set.txt"
            open(read_file, "wb").write(response.content)
            with open(read_file) as f:
                for line in f.readlines():
                    header = re.search(header_pattern, line)
                    if header:
                        read_entry = header.group(1)
                    content = re.search(content_pattern, read_entry)
                    if not content:
                        total_read += 1
                        if ((re.search(r'[infected|[iIu]mmune|challenge|tomi[zs]ed]',read_entry)) or (re.search(r'[Mm]is\w{0,2}RNA|lncRNA|circRNA|small RNA]',read_entry))):
                            unwanted_reads.append(read_entry)
                            continue
                        else:
                            wanted_reads.append(read_entry)
                            entry += 1
            if entry > 0:
                status = 'available'
        return status

    def get_assembly_metadata(self,id,gca):
        """ Function to retrieve assembly meta information from the public archives """
        try:
            handle = Entrez.esummary(db="assembly",id=id,report="full")
            record = Entrez.read(handle)
            meta_data = record['DocumentSummarySet']['DocumentSummary'][0]['Meta']
            bioproject_id = (record['DocumentSummarySet']['DocumentSummary'][0]['GB_BioProjects'][0]['BioprojectAccn']).strip() #This will return the bioproject id
            assembly_name = (record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyName']).strip()
            assembly_type = (record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyType']).strip()
            assembly_level = (record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyStatus']).strip()
            species_name = (record['DocumentSummarySet']['DocumentSummary'][0]['Organism']).strip()
            species_name = species_name.split('(')
            species_name = species_name[0].strip()
            assembly_path = (record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']).strip()
            assembly_path = assembly_path.replace('ftp:','https:');
            submitter = (record['DocumentSummarySet']['DocumentSummary'][0]['SubmitterOrganization']).strip()
            genome_rep = (record['DocumentSummarySet']['DocumentSummary'][0]['PartialGenomeRepresentation']).strip()
            taxonomy = (record['DocumentSummarySet']['DocumentSummary'][0]['Taxid']).strip()
            contig_N50 = (record['DocumentSummarySet']['DocumentSummary'][0]['ContigN50']).strip()
            scaffold_N50 = (record['DocumentSummarySet']['DocumentSummary'][0]['ScaffoldN50']).strip()
            contig_count = meta_data.split('<Stat category="contig_count" sequence_tag="all">')
            contig_count = contig_count[1].split('<')
            chromosome_count = meta_data.split('<Stat category="chromosome_count" sequence_tag="all">')
            chromosome_count = chromosome_count[1].split('<')
            scaffold_count = meta_data.split('<Stat category="scaffold_count" sequence_tag="all">')
            scaffold_count = scaffold_count[1].split('<')
            refseq_accession = (record['DocumentSummarySet']['DocumentSummary'][0]['Synonym']['RefSeq']).strip()
            assembly_date = (record['DocumentSummarySet']['DocumentSummary'][0]['SubmissionDate']).strip()
            genbank_accession = (record['DocumentSummarySet']['DocumentSummary'][0]['Synonym']['Genbank']).strip()
            assembly_date = assembly_date.split()
            meta_data = record['DocumentSummarySet']['DocumentSummary'][0]['Meta']
            total_length = meta_data.split('<Stat category="total_length" sequence_tag="all">')
            total_length = total_length[1].split('<')
            ungapped_length = meta_data.split('<Stat category="ungapped_length" sequence_tag="all">')
            ungapped_length = ungapped_length[1].split('<')
            total_gap_length = int(total_length[0]) - int(ungapped_length[0])
            wgs = (record['DocumentSummarySet']['DocumentSummary'][0]['WGS']).strip()
            if genome_rep == 'false':
                genome_rep = 'full'
            elif genome_rep == 'true':
                genome_rep = 'partial'  
        except HTTPError as e:
            content = e.read()
            raise eHive.CompleteEarlyException('Could not connect to download assembly meta data. Check that the assembly data exists in the public archives\n'+str(error))

        if assembly_level == 'Contig':
            top_level_count = contig_count[0]
        elif assembly_level == 'Scaffold':
            top_level_count = scaffold_count[0]
        elif assembly_level == 'Chromosome':
            top_level_count = chromosome_count[0]
        elif assembly_level == 'Complete Genome':
            top_level_count = scaffold_count[0]

        if genbank_accession == gca:
            found_record = 1
        else:
            found_record = 0
        #Summarise metadata information for reporting
        summary_report = [found_record,genbank_accession,wgs,bioproject_id,assembly_name,assembly_type,assembly_level,species_name,assembly_path,submitter,genome_rep,\
                          taxonomy,contig_N50,contig_count[0],chromosome_count[0],scaffold_count[0],scaffold_N50,refseq_accession,assembly_date[0],total_length[0],\
                          total_gap_length,top_level_count]
       
        return(summary_report)

    def set_assembly_group(self,bioproject):
        handle = Entrez.esearch(db="bioproject", term=bioproject)
        record = Entrez.read(handle)
        ids = record['IdList']
        ebp_dict = {}
        dtol_dict = {}
        vgp_dict = {}
        erga_dict = {}
        ergapp_dict = {}
        abp_dict = {}
        genomes25_dict = {}
        tol_dict = {}
        asg_dict = {}
        b10k_dict = {}
        g10k_dict = {}
        count = 0
        pri_asm_group = 'ungrouped' #default group status for new assemblies
        #Find all bioprojects under the Earth Genome  Project
        with os.popen("esearch -query 'PRJNA533106' -db bioproject | elink -related | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc ") as ebp:
            for id in ebp.readlines():
                count += 1
                ebp_dict[id] = count
        #Find all bioprojects under the G10K Project
        with os.popen("esearch -query 'PRJNA566188' -db bioproject | elink -related | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc ") as g10k:
            count = 0
            for id in g10k.readlines():
                count += 1
                g10k_dict[id] = count
        #Find all bioprojects under the Darwin Tree of Life Project
        with os.popen("esearch -query 'PRJEB43743' -db bioproject | elink -related | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc ") as asg:
            count = 0
            for id in asg.readlines():
                count += 1
                asg_dict[id] = count
        #Find all bioprojects under the Darwin Tree of Life Project
        with os.popen("esearch -query 'PRJEB40665' -db bioproject | elink -related | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc ") as dtol:
            count = 0
            for id in dtol.readlines():
                count += 1
                dtol_dict[id] = count
        #Find all bioprojects under the Tree of Life Project
        with os.popen("esearch -query 'PRJEB43745' -db bioproject | elink -related | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc ") as tol:
            count = 0
            for id in tol.readlines():
                count += 1
                tol_dict[id] = count
        #Find all bioprojects under the Vertebrates Genome Project
        with os.popen("esearch -query 'PRJNA489243' -db bioproject | elink -related | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc ") as vgp:
            count = 0
            for id in vgp.readlines():
                count += 1
                vgp_dict[id] = count
        #Find all bioprojects under the B10K  Project. Apparently, the B10K project is linked to two bioprojects
        with os.popen("esearch -query 'PRJNA489244' -db bioproject | elink -related | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc ") as b10k:
            count = 0
            for id in b10k.readlines():
                count += 1
                b10k_dict[id] = count
        with os.popen("esearch -query 'PRJNA545868' -db bioproject | elink -related | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc ") as b10:
            count = 0
            for id in b10.readlines():
                count += 1
                b10k_dict[id] = count
        #Find all bioprojects under the European Reference Genome Project
        with os.popen("esearch -query 'PRJEB43510' -db bioproject | elink -related | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc ") as erga:
            count = 0
            for id in erga.readlines():
                count += 1
                erga_dict[id] = count
        #Find all bioprojects under the African Biogenome Project
        with os.popen("esearch -query 'PRJNA811786' -db bioproject | elink -related | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc ") as abp:
            count = 0
            for id in abp.readlines():
                count += 1
                abp_dict[id] = count
        #Find all bioprojects under the European Reference Genome Atlas Pilot Project
        with os.popen("esearch -query 'PRJEB47820' -db bioproject | elink -related | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc ") as ergapp:
            count = 0
            for id in ergapp.readlines():
                count += 1
                ergapp_dict[id] = count
        #Find all bioprojects under the 25 Genomes Project
        with os.popen("esearch -query 'PRJEB33226' -db bioproject | elink -related | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc ") as gnomes25:
            count = 0
            for id in gnomes25.readlines():
                count += 1
                genomes25_dict[id] = count
        #Find linked projects for assembly being registered. Sometimes, assembly could be linked to multiple (related) projects
        with os.popen("esearch -query '"+bioproject.strip()+"' -db bioproject | elink -related | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc ") as bioproject:
            for id in bioproject.readlines():
                if id in b10k_dict.keys():
                    pri_asm_group = 'b10k'
                    return pri_asm_group
                elif id in dtol_dict.keys():
                    pri_asm_group = 'dtol'
                    return pri_asm_group
                elif id in erga_dict.keys():
                    pri_asm_group = 'erga'
                    return pri_asm_group
                elif id in vgp_dict.keys():
                    pri_asm_group = 'vgp'
                    return pri_asm_group
                elif id in abp_dict.keys():
                    pri_asm_group = 'abp'
                    return pri_asm_group
                elif id in ergapp_dict.keys():
                    pri_asm_group = 'ergapp'
                    return pri_asm_group
                elif id in genomes25_dict.keys():
                    pri_asm_group = 'genomes25'
                    return pri_asm_group
                elif id in tol_dict.keys():
                    pri_asm_group = 'tol'
                    return pri_asm_group
                elif id in asg_dict.keys():
                    pri_asm_group = 'asg'
                    return pri_asm_group
                elif id in ebp_dict.keys():
                    pri_asm_group = 'ebp'
                    return pri_asm_group
                else:
                    pri_asm_group = 'ungrouped'
                    return pri_asm_group
        return pri_asm_group
