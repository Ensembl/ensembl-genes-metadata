# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

""" Script to to check and resolve irregularities in meta database entries.  """
from datetime import datetime,timedelta,date
import pymysql
pymysql.install_as_MySQLdb()
import urllib.request,json
import argparse
import time
import sqlalchemy as db
import string
from pathlib import Path
import re
import errno
import os
from urllib.error import HTTPError
from Bio import Entrez
Entrez.email = os.getenv('genebuild_email')


def update_sec_asm_group_for_faang(db,host,port,user,password,group,taxon_id):
    #Update secondary assembly group for assemblies belonging to a particular species"
    con = pymysql.connect(
            host=host, user=user, passwd=password, port=int(port), db=db.strip(), charset='utf8mb4', cursorclass=pymysql.cursors.DictCursor
    )
    #get assembly id of entries matching species taxa
    sql = "SELECT assembly_id FROM assembly WHERE taxonomy = "+taxon_id
    info = fetch_data(sql,db,host,int(port),user,password)
    results = 0
    for row in info:
        try:
            with con.cursor() as cur:
                cur.execute('UPDATE meta SET sec_asm_group = %s WHERE assembly_id = %s' ,  (group, row))
        except Exception as error:
            print ('Oops! Secondary assembly group update did not complete successfully. Kindly check the error below\n'+str(error))
            break
        else:
            results += 1
            con.commit()
    print('Numer of assemblies with secondary group set to  '+group+' = '+str(results))

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

def get_ids(term):
    """ Check and return any linked assembly to accession of interest """
    ids = []
    handle = Entrez.esearch(db="assembly", term=term)
    record = Entrez.read(handle)
    ids.append(record["IdList"])
    return ids

def update_contigN50_and_scaffoldN50(db,host,port,user,password):
    con = pymysql.connect(
        host=host, user=user, passwd=password, port=int(port), db=db.strip(), charset='utf8mb4', cursorclass=pymysql.cursors.DictCursor
    )
    contig_N50_assemblies = {}
    scaffold_N50_assemblies = {}
    #Get assemblies with contig_N50 = 0 or scaffold_N50 = 0 for update
    contig_N50_query = "SELECT CONCAT(chain,'.',version),meta.assembly_id FROM assembly JOIN meta using(assembly_id) WHERE is_current = 1 AND contig_N50 = 0"
    scaffold_N50_query = "SELECT CONCAT(chain,'.',version),assembly_id FROM assembly JOIN meta using(assembly_id) WHERE is_current = 1 AND scaffold_N50 = 0 AND assembly_level IN ('scaffold','chromosome')"
    contig_N50_info = fetch_data(contig_N50_query,db,host,int(port),user,password)
    for row in contig_N50_info:
        contig_N50_assemblies[row[1]] = row[0]
    scaffold_N50_info = fetch_data(scaffold_N50_query,db,host,int(port),user,password)
    for row in scaffold_N50_info:
        scaffold_N50_assemblies[row[1]] = row[0]
    
    #Create a new dictionary that combines both contig_N50 and scaffold_N50 assemblies for update
    if len(contig_N50_assemblies) > 0:
        assemblies_to_update = contig_N50_assemblies.copy()
        assemblies_to_update.update(scaffold_N50_assemblies)
    elif len(scaffold_N50_assemblies) > 0:
        assemblies_to_update = scaffold_N50_assemblies.copy()
    else:
        #Terminate script as no assembly was found for update
        print('No assembly found for contig/scaffold update at this time')
        raise SystemExit()

    update_cnt = 0
    for key,val in assemblies_to_update.items():
        #Get contig_N50 and scaffold_N50 records
        contig_N50 = 0
        scaffold_N50 = 0
        for id in get_ids(val):
            if len(id) >= 1: #Assembly is linked to other assembly(ies)
                for i in range(len(id)): #Get metadata for each linked assembly
                    try:
                        handle = Entrez.esummary(db="assembly",id=id[i],report="full")
                        record = Entrez.read(handle, validate=False)
                        genbank_accession = (record['DocumentSummarySet']['DocumentSummary'][0]['Synonym']['Genbank']).strip()
                        meta_data = record['DocumentSummarySet']['DocumentSummary'][0]['Meta']
                        contig_N50 = (record['DocumentSummarySet']['DocumentSummary'][0]['ContigN50']).strip()
                        scaffold_N50 = (record['DocumentSummarySet']['DocumentSummary'][0]['ScaffoldN50']).strip()
                        if genbank_accession == val:
                            break
                    except Exception as error:
                        print ('Oops! Error fetching assembly report. Kindly check the error below\n'+str(error))
        try:
            #Check if any assembly needs update to contig/scaffold related fields
            if int(scaffold_N50) > 0:
                with con.cursor() as cur:
                    cur.execute('UPDATE meta SET contig_N50 = %s WHERE assembly_id = %s' ,  (contig_N50, key))
                    cur.execute('UPDATE meta SET scaffold_N50 = %s WHERE assembly_id = %s' ,  (scaffold_N50, key))
            elif int(contig_N50) > 0:
                #Update only contig related fields
                with con.cursor() as cur:
                    cur.execute('UPDATE meta SET contig_N50 = %s WHERE assembly_id = %s' ,  (contig_N50, key))
        except Exception as error:
            print ('Oops! Contig/Scaffold update did not complete successfully. Kindly check the error below\n'+str(error))
            break
            raise SystemExit()
        else:
            update_cnt += 1
            con.commit()
    print('Total number of assemblies updated = '+str(update_cnt))
    return

def get_rapid_dbs(release):
    #Get all dbnames from the vertebrates and metazoa divisions that are linked to the latest rapid release version
    rapid_dbs = {}
    db_meta = {}
    gmethod = ''
    accession = ''
    ini_date = ''
    last_date = ''
    count = 0

    # Processing rapid dbs
    sql = "SELECT gd.dbname FROM genome_database gd JOIN genome g ON g.genome_id = gd.genome_id JOIN organism o ON g.organism_id = o.organism_id JOIN data_release dr ON g.data_release_id = dr.data_release_id JOIN division d ON g.division_id = d.division_id WHERE dr.ensembl_genomes_version = " + release + " AND d.name IN ('EnsemblVertebrates','EnsemblMetazoa') AND gd.type = 'core' ORDER by gd.dbname"
    results = fetch_data(sql,os.getenv('ENS_META_QRP'),os.getenv('PROD_META'),int(os.getenv('PROD_PORT')),os.getenv('GBUSER_R'),'')
    for row in results:
        rapid_dbs[row[0]] = release
        count+=1

    # Queries to fetch needed meta table info 
    assembly_accession_query = "SELECT meta_value FROM meta WHERE meta_key = 'assembly.accession'"
    initial_rel_date_query = "SELECT meta_value FROM meta WHERE meta_key = 'genebuild.initial_release_date'"
    last_geneset_update_query = "SELECT meta_value FROM meta WHERE meta_key = 'genebuild.last_geneset_update'"
    end_date = ''

    #Retrieve metadata from each rapid db to use in updating the meta database table (genebuild_status)
    for row in rapid_dbs:
        if (re.search(r'lingula_anatina_gca001039355v2_core',str(row))): #This db entry exists in the production meta database but it is not live
            continue
        # Rapid live servers is queried based on the release version
        if int(release) % 2 == 0:
            #Version is even, so query liver server 2
            assembly_accession = fetch_data(assembly_accession_query,row,os.getenv('RRHOST2'),int(os.getenv('RRWPORT')),os.getenv('GBUSER_R'),'')
            for acc in assembly_accession:
                accession = acc[0]
            initial_rel_date = fetch_data(initial_rel_date_query,row,os.getenv('RRHOST2'),int(os.getenv('RRWPORT')),os.getenv('GBUSER_R'),'')
            # End date can be fetched in either ways depending on what is set at the time
            if initial_rel_date:
                for ind in initial_rel_date:
                    last_date = ind[0]
            else:
                last_geneset_update = fetch_data(last_geneset_update_query,row,os.getenv('RRHOST2'),int(os.getenv('RRWPORT')),os.getenv('GBUSER_R'),'')
                for ld in last_geneset_update:
                    last_date = ld[0]
            genebuild_start_date_query = "SELECT date_started FROM genebuild_status WHERE assembly_accession = '"+accession+"' AND is_current = 1"
            genebuild_start_date = fetch_data(genebuild_start_date_query,os.getenv('REG_DB'),os.getenv('GBS1'),int(os.getenv('GBP1')),os.getenv('GBUSER_R'),'')
            for gsd in genebuild_start_date:
                start_date = gsd[0]
        else:
            # Query liver server 1
            assembly_accession = fetch_data(assembly_accession_query,row,os.getenv('RRHOST1'),int(os.getenv('RRWPORT')),os.getenv('GBUSER_R'),'')
            for acc in assembly_accession:
                accession = acc[0]
            initial_rel_date = fetch_data(initial_rel_date_query,row,os.getenv('RRHOST1'),int(os.getenv('RRWPORT')),os.getenv('GBUSER_R'),'')
            # End date can be fetched in either ways depending on what is set at the time
            if initial_rel_date:
                for ind in initial_rel_date:
                    last_date = ind[0]
            else:
                last_geneset_update = fetch_data(last_geneset_update_query,row,os.getenv('RRHOST1'),int(os.getenv('RRWPORT')),os.getenv('GBUSER_R'),'')
                for ld in last_geneset_update:
                    last_date = ld[0]
            genebuild_start_date_query = "SELECT date_started FROM genebuild_status WHERE assembly_accession = '"+accession+"' AND is_current = 1"
            genebuild_start_date = fetch_data(genebuild_start_date_query,os.getenv('REG_DB'),os.getenv('GBS1'),int(os.getenv('GBP1')),os.getenv('GBUSER_R'),'')
            for gsd in genebuild_start_date:
                start_date = gsd[0]

        # Format the genebuild_start_date and end_dates to match database type
        # This is an estimate of 3 weeks since we dont actually store the day part in the meta table
        #start_date = re.sub('EnsemblMetazoa|Ensembl|WormBase|Plants', '01', start_date)
        if last_date:
            last_date = last_date + '-21'
        else:
                
            last_date = start_date + timedelta(days = 21)

        db_meta[accession] = accession + '\t' + str(start_date) + '\t' + str(last_date) + '\t' + row
        last_date = ''
    print('Rapid dbs = '+str(len(db_meta)))
    return db_meta

""" Method to fetch current live dbs from the Main Release Site """
def get_main_dbs(release):
    #Get all dbnames from the vertebrates and metazoa divisions that are linked to the latest main release version
    main_dbs = {}
    db_meta = {}
    accession = ''
    ini_date = ''
    last_date = ''
    start_date = ''
    count = 0

    # Processing main dbs
    sql = "SELECT gd.dbname FROM genome_database gd JOIN genome g ON g.genome_id = gd.genome_id JOIN organism o ON g.organism_id = o.organism_id JOIN data_release dr ON g.data_release_id = dr.data_release_id JOIN division d ON g.division_id = d.division_id WHERE dr.ensembl_version = " + release + " AND d.name IN ('EnsemblVertebrates') AND gd.type = 'core' ORDER by gd.dbname"
    results = fetch_data(sql,os.getenv('ENS_META'),os.getenv('PROD_META'),int(os.getenv('PROD_PORT')),os.getenv('GBUSER_R'),'')
    for row in results:
        main_dbs[row[0]] = release
        count+=1
    # Queries to fetch needed meta table info
    assembly_accession_query = "SELECT meta_value FROM meta WHERE meta_key = 'assembly.accession'"
    initial_rel_date_query = "SELECT meta_value FROM meta WHERE meta_key = 'genebuild.initial_release_date'"
    last_geneset_update_query = "SELECT meta_value FROM meta WHERE meta_key = 'genebuild.last_geneset_update'"
    end_date = ''

    #Retrieve metadata from each main db to use in updating the meta database table (genebuild_status)
    for row in main_dbs:
        assembly_accession = fetch_data(assembly_accession_query,row,os.getenv('MIRROR'),int(os.getenv('MPORT')),os.getenv('GBUSER_R'),'')
        for acc in assembly_accession:
            accession = acc[0]
        initial_rel_date = fetch_data(initial_rel_date_query,row,os.getenv('MIRROR'),int(os.getenv('MPORT')),os.getenv('GBUSER_R'),'')
        # End date can be fetched in either ways depending on what is set at the time
        if initial_rel_date:
            for ind in initial_rel_date:
                last_date = ind[0]
        else:
            last_geneset_update = fetch_data(last_geneset_update_query,row,os.getenv('MIRROR'),int(os.getenv('MPORT')),os.getenv('GBUSER_R'),'')
            for ld in last_geneset_update:
                last_date = ld[0]
        genebuild_start_date_query = "SELECT date_started FROM genebuild_status WHERE assembly_accession = '"+accession+"' AND is_current = 1"
        genebuild_start_date = fetch_data(genebuild_start_date_query,os.getenv('REG_DB'),os.getenv('GBS1'),int(os.getenv('GBP1')),os.getenv('GBUSER_R'),'')
        for gsd in genebuild_start_date:
            start_date = gsd[0]

        # Format the genebuild_start_date and end_dates to match database type
        # This is an estimate of 3 weeks since we dont actually store the day part in the meta table
        if last_date:
            last_date = last_date + '-21'
        else:

            last_date = start_date + timedelta(days = 21)
        db_meta[accession] = accession +'\t' + str(start_date) + '\t' + str(last_date) + '\t' + row
        last_date = ''
    print('Main dbs = '+str(len(db_meta)))
    return db_meta

def update_genebuild_status(db,host,port,user,password):
    #Function to update genebuild status columns where db has been released but related fields don't reflect the full information
    con = pymysql.connect(
        host=host, user=user, passwd=password, port=int(port), db=db.strip(), charset='utf8mb4', cursorclass=pymysql.cursors.DictCursor
    )
    rapid_meta = {}
    main_meta = {}
    version = ''
    release = get_current_release()
    if release:
        version = release.split('\t')
        print('Live rapid release version = '+str(version[2])+ ' and release date is '+str(version[3]))
        print('Main release version = '+str(version[0])+ ' and release date is '+str(version[1]))
        rapid_meta = get_rapid_dbs(version[2])
        main_meta = get_main_dbs(version[0])

    # Merge both rapid and main dbs for comparison
    live_db_meta_info = main_meta.copy()
    live_db_meta_info.update(rapid_meta)

    #Get assemblies that have been released but with incorrect or incomplete genebuild status
    genebuild_status_query = "SELECT assembly_accession,date_completed,release_server,release_date,release_version FROM genebuild_status WHERE progress_status != 'handed over' OR release_server = 'not available' AND annotation_source IN ('ensembl','draft') AND is_current = 1"
    genebuild_status = {}
    results = fetch_data(genebuild_status_query,db,host,int(port),user,password)
    for row in results:
        genebuild_status[row[0]] = row
    print('Incomplete genebuild entries= '+str(len(genebuild_status)))
    date_cnt = 0
    prog_cnt = 0 
    cnnt = 0
    # Search all current genebuild entries to match against the meta info retrieved from the live servers
    for key,val in genebuild_status.items():
        if live_db_meta_info.get(key):
            line = live_db_meta_info[key]
            meta_info = line.split('\t')
            genebuild_meta = genebuild_status[meta_info[0]]
            # Does accession exist in registry?
            if not meta_info[0] in genebuild_status.keys():
                continue
            else:
                try:
                    with con.cursor() as cur:
                        # Check to ensure only annotations done by the genebuild team are updated
                        if genebuild_meta[1] is None:
                            date_cnt += 1
                            cur.execute('UPDATE genebuild_status set date_completed =  %s WHERE assembly_accession = %s AND is_current = %s',(meta_info[2], meta_info[0], 1))
                            cur.execute('UPDATE genebuild_status set progress_status =  %s WHERE assembly_accession = %s AND is_current = %s',('handed over', meta_info[0], 1))
                        else:
                            prog_cnt += 1
                            cur.execute('UPDATE genebuild_status set progress_status =  %s WHERE assembly_accession = %s AND is_current = %s',('handed over', meta_info[0], 1))

                            # Where assembly has been released on both rapid and main, set server to Main as rapid would soon fade away
                        if rapid_meta.get(key) and main_meta.get(key):
                            cur.execute('UPDATE genebuild_status set release_version =  %s WHERE assembly_accession = %s AND is_current = %s',(version[0],meta_info[0], 1))
                            cur.execute('UPDATE genebuild_status set release_date =  %s WHERE assembly_accession = %s AND is_current = %s',(version[1],meta_info[0], 1))
                            cur.execute('UPDATE genebuild_status set release_server =  %s WHERE assembly_accession = %s AND is_current = %s',('main',meta_info[0], 1))
                        elif not rapid_meta.get(meta_info[0]):
                            #Main only
                            cur.execute('UPDATE genebuild_status set release_version =  %s WHERE assembly_accession = %s AND is_current = %s',(version[0],meta_info[0], 1))
                            cur.execute('UPDATE genebuild_status set release_date =  %s WHERE assembly_accession = %s AND is_current = %s',(version[1],meta_info[0], 1))
                            cur.execute('UPDATE genebuild_status set release_server =  %s WHERE assembly_accession = %s AND is_current = %s',('main',meta_info[0], 1))
                        else:
                            #Rapid only
                            cur.execute('UPDATE genebuild_status set release_version =  %s WHERE assembly_accession = %s AND is_current = %s',(version[2],meta_info[0], 1))
                            cur.execute('UPDATE genebuild_status set release_date =  %s WHERE assembly_accession = %s AND is_current = %s',(version[3],meta_info[0], 1))
                            cur.execute('UPDATE genebuild_status set release_server =  %s WHERE assembly_accession = %s AND is_current = %s',('rapid',meta_info[0], 1))

                except Exception as error:
                    print ('Genebuild update did not complete successfully. Kindly check the error below\n'+str(error))
                    break
                    raise SystemExit()
                else:
                    con.commit()
        else:
            cnnt += 1
            with open('rpt.txt', 'a') as writer:
                writer.write(str(val))
            writer.close()
    print('Assembly not in live severs '+str(cnnt))
    print('Genebuild updated with completion date = '+ str(date_cnt))
    print('Genebuild progress status update = '+ str(prog_cnt))
    return

""" This method retrieves the current Ensembl release version for both Main and Rapid releases """
def get_current_release():
    main_release = 0
    main_rel_date = ''
    #Obtain release number for main
    sql = "SELECT ensembl_version,release_date FROM data_release WHERE is_current = 1 ORDER BY ensembl_version DESC LIMIT 1"
    results = fetch_data(sql,os.getenv('ENS_META'),os.getenv('PROD_META'),int(os.getenv('PROD_PORT')),os.getenv('GBUSER_R'),'')
    for row in results:
        # In order to find the cuorrect set of db on liver servers, we use one version less the current release
        main_release = row[0] - 1
        main_rel_date = row[1]

    #Obtain release number for RR
    rr_release = 0
    rr_rel_date = ''
    sql = "SELECT ensembl_genomes_version,release_date FROM data_release WHERE is_current = 1 ORDER BY ensembl_genomes_version DESC LIMIT 1"
    results = fetch_data(sql,os.getenv('ENS_META_QRP'),os.getenv('PROD_META'),int(os.getenv('PROD_PORT')),os.getenv('GBUSER_R'),'')
    for row in results:
        # In order to find the cuorrect set of db on liver servers, we use one version less the current release
        rr_release = row[0] - 1
        rr_rel_date = row[1]
    results = str(main_release) + '\t' + str(main_rel_date) + '\t' + str(rr_release) + '\t' + str(rr_rel_date)

    return results

def update_pri_assembly_group(db,host,port,user,password,accession_list):
    #Function update incorrectly assigned assembly group
    con = pymysql.connect(
            host=host, user=user, passwd=password, port=int(port), db=db.strip(), charset='utf8mb4', cursorclass=pymysql.cursors.DictCursor
    )
    bioproject_id = ''
    with open(cfg_path) as file:
        while True:
            entry = file.readline()
            if not entry:
                break
            accession = entry.split('\t')
            accession[0].replace('\n','').strip()
            print('Accession being processed is '+str(accession[0].strip()))
            with os.popen("esearch -query "+accession[0].strip()+" -db assembly | elink -target bioproject | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc ") as bioproject:
                for bioproject_id in bioproject.readlines():
                    bioproject_id.strip()

            results = 0
            ebp_dict = {}
            ebp_cnt = 0
            argp_dict = {}
            argp_cnt = 0
            dtol_dict = {}
            dtol_cnt = 0
            vgp_dict = {}
            vgp_cnt = 0
            erga_dict = {}
            erga_cnt = 0
            ergapp_dict = {}
            ergapp_cnt = 0
            abp_dict = {}
            abp_cnt = 0
            genomes25_dict = {}
            genomes25_cnt = 0
            tol_dict = {}
            tol_cnt = 0
            asg_dict = {}
            asg_cnt = 0
            b10k_dict = {}
            b10k_cnt = 0
            g10k_dict = {}
            g10k_cnt = 0
            cbp_dict = {}
            cbp_cnt = 0
            count = 0
            ungrouped_cnt = 0
            #Secondary assembly group settings
            asm_group_dict = {'cbp' : 'ebp', 'argp' : 'ebp', 'genomes25' : 'ebp', 'ebp' : 'None', 'abp' : 'ebp', 'asg' : 'ebp', 'b10k' : 'ebp', 'dtol' : 'erga ebp', 'erga' : 'ebp', 'ergapp' : 'erga', 'g10k' : 'erga  ebp', 'tol': 'None', 'ungrouped' : 'None', 'vgp' : 'ebp'}

            pri_asm_group = 'ungrouped' #default group status for new assemblies
            sec_asm_group = 'None'
            #Find all bioprojects under the Canadian BioGenome  Project
            with os.popen("esearch -query 'PRJNA813333' -db bioproject | elink -related | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc ") as cbp:
                count = 0
                for id in cbp.readlines():
                    count += 1
                    cbp_dict[id.strip()] = count
            #Find all bioprojects under the Earth Genome  Project
            with os.popen("esearch -query 'PRJNA533106' -db bioproject | elink -related | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc ") as ebp:
                count = 0
                for id in ebp.readlines():
                    count += 1
                    ebp_dict[id.strip()] = count
                #Find all bioprojects under the Anopheles Genome  Project
            with os.popen("esearch -query 'PRJEB51690' -db bioproject | elink -related | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc ") as argp:
                count = 0
                for id in argp.readlines():
                    count += 1
                    argp_dict[id.strip()] = count
            #Find all bioprojects under the G10K Project
            with os.popen("esearch -query 'PRJNA566188' -db bioproject | elink -related | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc ") as g10k:
                count = 0
                for id in g10k.readlines():
                    count += 1
                    g10k_dict[id.strip()] = count
            #Find all bioprojects under the Darwin Tree of Life Project
            with os.popen("esearch -query 'PRJEB43743' -db bioproject | elink -related | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc ") as asg:
                count = 0
                for id in asg.readlines():
                    count += 1
                    asg_dict[id.strip()] = count
            #Find all bioprojects under the Darwin Tree of Life Project
            with os.popen("esearch -query 'PRJEB40665' -db bioproject | elink -related | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc ") as dtol:
                count = 0
                for id in dtol.readlines():
                    count += 1
                    dtol_dict[id.strip()] = count
            #Find all bioprojects under the Tree of Life Project
            with os.popen("esearch -query 'PRJEB43745' -db bioproject | elink -related | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc ") as tol:
                count = 0
                for id in tol.readlines():
                    count += 1
                    tol_dict[id.strip()] = count
            #Find all bioprojects under the Vertebrates Genome Project
            with os.popen("esearch -query 'PRJNA489243' -db bioproject | elink -related | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc ") as vgp:
                count = 0
                for id in vgp.readlines():
                    count += 1
                    vgp_dict[id.strip()] = count
            #Find all bioprojects under the B10K  Project. Apparently, the B10K project is linked to two bioprojects
            with os.popen("esearch -query 'PRJNA489244' -db bioproject | elink -related | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc ") as b10k:
                count = 0
                for id in b10k.readlines():
                    count += 1
                    b10k_dict[id.strip()] = count
            with os.popen("esearch -query 'PRJNA545868' -db bioproject | elink -related | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc ") as b10:
                count = 0
                for id in b10.readlines():
                    count += 1
                    b10k_dict[id.strip()] = count
            #Find all bioprojects under the European Reference Genome Project
            with os.popen("esearch -query 'PRJEB43510' -db bioproject | elink -related | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc ") as erga:
                count = 0
                for id in erga.readlines():
                    count += 1
                    erga_dict[id.strip()] = count
            #Find all bioprojects under the African Biogenome Project
            with os.popen("esearch -query 'PRJNA811786' -db bioproject | elink -related | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc ") as abp:
                count = 0
                for id in abp.readlines():
                    count += 1
                    abp_dict[id.strip()] = count
            #Find all bioprojects under the European Reference Genome Atlas Pilot Project
            with os.popen("esearch -query 'PRJEB47820' -db bioproject | elink -related | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc ") as ergapp:
                count = 0
                for id in ergapp.readlines():
                    count += 1
                    ergapp_dict[id.strip()] = count
            #Find all bioprojects under the 25 Genomes Project
            with os.popen("esearch -query 'PRJEB33226' -db bioproject | elink -related | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc ") as gnomes25:
                count = 0
                for id in gnomes25.readlines():
                    count += 1
                    genomes25_dict[id.strip()] = count
            #Find linked projects for assembly being registered. Sometimes, assembly could be linked to multiple (related) projects
            with os.popen("esearch -query '"+bioproject_id.strip()+"' -db bioproject | elink -related | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc ") as bioproject:
                for id in bioproject.readlines():
                    
                    if id.strip() in dtol_dict.keys() and id.strip() in vgp_dict.keys():
                        #Where assembly is linked to VGP and DToL  projects, prioritise DToL over VGP
                        pri_asm_group = 'dtol'
                        dtol_cnt += 1
                        break
                    elif id.strip() in dtol_dict.keys() and id.strip() in erga_dict.keys():
                        #Where assembly is linked to ERGA and DToL  projects, prioritise DToL over ERGA
                        pri_asm_group = 'dtol'
                        dtol_cnt += 1
                        break
                    elif id.strip() in dtol_dict.keys() and id.strip() in ebp_dict.keys():
                        #Where assembly is linked to EBP and DToL  projects, prioritise DToL over EBP
                        pri_asm_group = 'dtol'
                        dtol_cnt += 1
                        break
                    elif id.strip() in vgp_dict.keys() and id.strip() in ergapp_dict.keys():
                        #Where assembly is linked to ERGAPP and VGP, prioritise ERGA over others
                        pri_asm_group = 'erga'
                        erga_cnt += 1
                        break
                    elif id.strip() in vgp_dict.keys() and id.strip() in erga_dict.keys():
                        #Where assembly is linked to ERGA and VGP, prioritise ERGA over VGP
                        pri_asm_group = 'erga'
                        erga_cnt += 1
                        break
                    elif id.strip() in erga_dict.keys() and id.strip() in ebp_dict.keys():
                        #Where assembly is linked to ERGA and EBP, prioritise ERGA over EBP
                        pri_asm_group = 'erga'
                        erga_cnt += 1
                        break
                    elif id.strip() in vgp_dict.keys() and id.strip() in ebp_dict.keys():
                        #Where assembly is linked to VGP and its related projects, prioritise VGP over others
                        pri_asm_group = 'vgp'
                        vgp_cnt += 1
                        break
                    elif id.strip() in vgp_dict.keys() and id.strip() in g10k_dict.keys():
                        #Where assembly is linked to VGP and its related project G10K, prioritise VGP over G10K
                        pri_asm_group = 'vgp'
                        vgp_cnt += 1
                        break
                    elif id.strip() in vgp_dict.keys() and id.strip() in b10k_dict.keys():
                        #Where assembly is linked to VGP and its related projects, prioritise VGP over others
                        pri_asm_group = 'vgp'
                        vgp_cnt += 1
                        break
                    elif id.strip() in vgp_dict.keys():
                        pri_asm_group = 'vgp'
                        vgp_cnt += 1
                        break
                    elif id.strip() in abp_dict.keys():
                        pri_asm_group = 'abp'
                        abp_cnt += 1
                        break
                    elif id.strip() in argp_dict.keys():
                        pri_asm_group = 'argp'
                        argp_cnt += 1
                        break
                    elif id.strip() in b10k_dict.keys():
                        pri_asm_group = 'b10k' 
                        b10k_cnt += 1
                        break
                    elif id.strip() in g10k_dict.keys():
                        pri_asm_group = 'g10k'
                        g10k_cnt += 1
                        break
                    elif id.strip() in cbp_dict.keys():
                        pri_asm_group = 'cbp'
                        cbp_cnt += 1
                        break
                    elif id.strip() in dtol_dict.keys():
                        pri_asm_group = 'dtol'
                        dtol_cnt += 1
                        break
                    elif id.strip() in ergapp_dict.keys():
                        pri_asm_group = 'ergapp'
                        ergapp_cnt += 1
                        break
                    elif id.strip() in erga_dict.keys():
                        pri_asm_group = 'erga'
                        erga_cnt += 1
                        break
                    elif id.strip() in genomes25_dict.keys():
                        pri_asm_group = 'genomes25'
                        genomes25_cnt += 1
                        break
                    elif id.strip() in tol_dict.keys():
                        pri_asm_group = 'tol'
                        tol_cnt += 1
                        break
                    elif id.strip() in asg_dict.keys():
                        pri_asm_group = 'asg'
                        asg_cnt += 1
                        break
                    elif id.strip() in ebp_dict.keys():
                        pri_asm_group = 'ebp'
                        ebp_cnt += 1
                        break
                    else:
                        pri_asm_group = 'ungrouped'
                        ungrouped_cnt += 1
                      #  break
            if asm_group_dict[pri_asm_group]:
                sec_asm_group = asm_group_dict[pri_asm_group]
            try:
                with con.cursor() as cur:
                    #print('Attempt to update db with group = '+bioproject_id +'\t'+id+'\t'+pri_asm_group)
                    cur.execute('UPDATE assembly SET pri_asm_group = %s WHERE CONCAT(chain,".",version) = %s' ,  (pri_asm_group, accession[0]))
                    cur.execute('UPDATE meta SET sec_asm_group = %s WHERE assembly_id = %s' ,  (sec_asm_group, accession[1]))
            except Exception as error:
                print ('Oops! Primary assembly group update did not complete successfully. Kindly check the error below\n'+str(error))
                break
            else:
                con.commit()
            print('Aassembly group assigned to accession '+accession[0]+' is Primary: '+pri_asm_group+' Secondary: '+sec_asm_group)
            
    report = 'DToL = '+str(dtol_cnt)+'\nVGP = '+str(vgp_cnt)+'\nERGA = '+str(erga_cnt)+'\nCBP = '+str(cbp_cnt)+'\nEBP = '+str(ebp_cnt)+'\nARGP = '+str(argp_cnt)+'\nASG = '+str(asg_cnt)+'\nB10K = '+str(b10k_cnt)+'\nG10K = '+str(g10k_cnt)+'\nGenomes 25 = '+str(genomes25_cnt)+'\nToL = '+str(tol_cnt)+'\nERGAP Pilot = '+str(ergapp_cnt)+'\nABP = '+str(abp_cnt)+'\n'
    print('Report of assembly group update run = '+str(report))
    print('Number of assemblies still classed as ungrouped = '+str(ungrouped_cnt))

    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--port', help='Port number for host', required=True)
    parser.add_argument('--password', help='Password for host', required=False)
    parser.add_argument('--user', help='Mysql user', required=True)
    parser.add_argument('--dbname', help='Database to be backed up', required=True)
    parser.add_argument('--server', help='Host server for database', required=True)
    parser.add_argument('--accessions', help='List of accessions to process', required=False)
    parser.add_argument('--asm_group', help='Assembly group', required=False)
    args = parser.parse_args()
    host = args.server
    db = args.dbname
    user = args.user
    port = args.port
    accession_list = args.accessions
    group = args.asm_group
    password = args.password
    approved_tasks = ['1 - update_sec_asm_group_for_faang','2 - update_pri_assembly_group','3 - update_contigN50_and_scaffoldN50','4 - update_genebuild_status']

    print('Enter the number corresponding to the task you wish to perfom and hit the Enter key\n')
    for task in approved_tasks:
        print(task+'\n')
    try:
        task = input('Enter the number matching the task you wish to perform above\n')
    except Exception as EOFError:
        raise ValueError('No response was received. Registration terminated.\n')
    if ((re.search(r'update_sec_asm_group_for_faang',str(task))) or task == '1'):
        update_sec_asm_group_for_faang(db,host,port,user,password,group,taxon_id)
    elif ((re.search(r'update_pri_assembly_group',str(task))) or task == '2'):
        update_pri_assembly_group(db,host,port,user,password,accession_list)
    elif ((re.search(r'update_contigN50_and_scaffoldN50',str(task))) or task == '3'):
        update_contigN50_and_scaffoldN50(db,host,port,user,password)
    elif ((re.search(r'update_genebuild_status',str(task))) or task == '6'):
        update_genebuild_status(db,host,port,user,password)
    else:
        print('No valid selection was made, hence script termination')

