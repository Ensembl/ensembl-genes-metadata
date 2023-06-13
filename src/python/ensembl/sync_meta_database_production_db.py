""" This script will query the production meta database to find all Vertebrates and Metazoa division dbs associated with the current/most recent release
    The  meta table of each db returned is queried to obtain key annotation parameters to update the genebuild status  of the assembly in the meta database
"""

import re
import os
import urllib.request,json
from urllib.error import HTTPError
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import text
import sys
from Bio import Entrez
import requests
Entrez.email = os.getenv('genebuild_email')

import pymysql
pymysql.install_as_MySQLdb()

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

""" This method takes a couple of parameters and  updates the annotation status of an assembly in the registry """
def update_genebuild_status(current_genebuild,live_db_meta_info,release_version,rel_date,host,database,user,port,password):
    con = pymysql.connect(
            host=host, user=user, passwd=password, port=port, database=database.strip(), charset='utf8mb4', cursorclass=pymysql.cursors.DictCursor
    )
    update_cnt = 0
    # Search all current genebuild entries to match against the meta info retrieved from the live servers
    for key,val in current_genebuild.items():
        if live_db_meta_info.get(key):
            #for row in live_db_meta_info:
            line = live_db_meta_info[key]
            meta_info = line.split('\t')
            # Does accession exist in registry?
            if not meta_info[0] in current_genebuild.keys():
                continue
            else:
                genebuild_meta = current_genebuild[meta_info[0]]
                if genebuild_meta[2] != 'handed over':
                    status = 'handed over'
                try:
                    with con.cursor() as cur:
                        update_report = meta_info[0]+'\t'+genebuild_meta+'\t'+meta_info[1]+'\n'
                        # Check to ensure only annotations done by the genebuild team are updated
                        if (genebuild_meta != meta_info[1] and meta_info[1].find('import') < 0):
                            update_cnt += 1
                            cur.execute('UPDATE genebuild_status set annotation_source =  %s WHERE assembly_accession = %s AND is_current = %s',(meta_info[1], meta_info[0], 1))
                            print('Progress status is '+status+' and accession is '+meta_info[0])
                            cur.execute('UPDATE genebuild_status set progress_status =  %s WHERE assembly_accession = %s AND is_current = %s',(status, meta_info[0], 1))
                            # At the moment, we are storing release information into the release_history table as well.
                            # This is temporary and should be removed once the team decides the table is no longer useful 
                            cur.execute('INSERT INTO release_history VALUES (%s,%s,%s,%s,%s)',(release_version,meta_info[0],meta_info[4],'rapid',rel_date))
                            cur.execute('UPDATE genebuild_status set release_version =  %s WHERE assembly_accession = %s AND is_current = %s',(release_version,meta_info[0], 1))
                            cur.execute('UPDATE genebuild_status set release_date =  %s WHERE assembly_accession = %s AND is_current = %s',(rel_date,meta_info[0], 1))
                            cur.execute('UPDATE genebuild_status set release_server =  %s WHERE assembly_accession = %s AND is_current = %s',('rapid',meta_info[0], 1))
                            # Update completion date if missing/incorrectly set
                            if meta_info[2] is None:
                                cur.execute('UPDATE genebuild_status set date_completed =  %s WHERE assembly_accession = %s AND is_current = %s',(meta_info[3], meta_info[0], 1))
                            with open('genebuild_update', 'a') as writer:
                                writer.write(update_report+meta_info[0]+'\t'+meta_info[4]+'\trapid\t'+str(rel_date)+'\t'+release_version+'\n')
                            writer.close()
                        else:
                            with open('no_update', 'a') as writer:
                                writer.write(update_report)
                            writer.close()
                            continue
        
                except Exception as error:
                    print(meta_info[0]+'\tCurrent status: '+genebuild_meta+'\tLive status: '+meta_info[1])
                    print ('Oops! Genebuild update did not complete successfully. Kindly check the error below\n'+str(error))
                    break
                else:
                    con.commit()
    print('Genebuild updated = '+ str(update_cnt))
        

""" This method retrieves all genebuild entries from the database """
def get_pending_genebuilds_for_update():
    #Get all current genebuild entries from the meta database that can be updated
    existing_accessions = {}
    sql = "SELECT assembly_accession, annotation_source, date_completed FROM genebuild_status WHERE is_current = 1 AND annotation_source = 'pending'"
    results = fetch_data(sql,os.getenv('REG_DB'),os.getenv('GBS1'),int(os.getenv('GBP1')),os.getenv('GBUSER_R'),'')
    for row in results:
        existing_accessions[row[0]] = row[1]
    print('Pending genebuilds = '+str(len(existing_accessions)))
    return existing_accessions

""" Method to fetch current live dbs from the Rapid Release Site """
def get_rapid_dbs(release):
    #Get all dbnames from the vertebrates and metazoa divisions that are linked to the latest rapid release version
    dbs = {}
    db_meta = {}
    gmethod = ''
    accession = ''
    prefix = ''
    ini_date = ''
    last_date = ''
    start_date = ''
    #release = release.split('\t')
    count = 0

    sql = "SELECT gd.dbname FROM genome_database gd JOIN genome g ON g.genome_id = gd.genome_id JOIN organism o ON g.organism_id = o.organism_id JOIN data_release dr ON g.data_release_id = dr.data_release_id JOIN division d ON g.division_id = d.division_id WHERE dr.ensembl_genomes_version = " + release + " AND d.name IN ('EnsemblVertebrates','EnsemblMetazoa') AND gd.type = 'core' ORDER by gd.dbname"
    results = fetch_data(sql,os.getenv('ENS_META_QRP'),os.getenv('PROD_META'),int(os.getenv('PROD_PORT')),os.getenv('GBUSER_R'),'')
    for row in results:
        dbs[row[0]] = release
        count+=1
    # Queries to fetch needed meta table info
    genebuild_method_query = "SELECT meta_value FROM meta WHERE meta_key = 'genebuild.method'"
    assembly_accession_query = "SELECT meta_value FROM meta WHERE meta_key = 'assembly.accession'"
    species_prefix_query = "SELECT meta_value FROM meta WHERE meta_key = 'species.stable_id_prefix'"
    initial_rel_date_query = "SELECT meta_value FROM meta WHERE meta_key = 'genebuild.initial_release_date'"
    last_geneset_update_query = "SELECT meta_value FROM meta WHERE meta_key = 'genebuild.last_geneset_update'"
    genebuild_start_date_query = "SELECT meta_value FROM meta WHERE meta_key = 'genebuild.start_date'"
    end_date = ''

    #Retrieve metadata from each db to use in updating the meta database table (genebuild_status)
    for row in dbs:
        if (re.search(r'lingula_anatina_gca001039355v2_core',str(row))): #This db entry exists in the production meta database but it is not live
            continue
        # Rapid live servers is queried based on the release version
        if int(release[1]) % 2 == 0:
            #Version is even, so query liver server 2
            genebuild_method = fetch_data(genebuild_method_query,row,os.getenv('RRHOST2'),int(os.getenv('RRWPORT')),os.getenv('GBUSER_R'),'')
            for mtd in genebuild_method:
                gmethod = mtd[0]
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
            genebuild_start_date = fetch_data(genebuild_start_date_query,row,os.getenv('RRHOST2'),int(os.getenv('RRWPORT')),os.getenv('GBUSER_R'),'')
            for gsd in genebuild_start_date:
                start_date = gsd[0]
        else:
            # Query liver server 1
            genebuild_method = fetch_data(genebuild_method_query,row,os.getenv('RRHOST1'),int(os.getenv('RRWPORT')),os.getenv('GBUSER_R'),'')
            for mtd in genebuild_method:
                gmethod = mtd[0] 
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
            genebuild_start_date = fetch_data(genebuild_start_date_query,row,os.getenv('RRHOST1'),int(os.getenv('RRWPORT')),os.getenv('GBUSER_R'),'')
            for gsd in genebuild_start_date:
                start_date = gsd[0]

            # Format the genebuild_start_date and end_dates to match database type
            # This is an estimate of 3 weeks since we dont actually store the day part in the meta table
        start_date = re.sub('EnsemblMetazoa|Ensembl|WormBase|Plants', '01', start_date)
        if last_date:
            last_date = last_date + '-21'
        else:
            last_date = re.sub('-01', '-21', start_date)

        # Store retrieved meta table info for each database
        main_methods = ['full_genebuild', 'anno']
        alt_methods = ['projection_build', 'braker']
        if any(c in gmethod for c in main_methods):
            gmethod = 'ensembl'
        elif any(c in gmethod for c in alt_methods):
            gmethod = 'draft'
        db_meta[accession] = accession + '\t' + gmethod + '\t' + start_date + '\t' + last_date + '\t' + row
        gmethod = ''
        last_date = ''
    print('Live dbs = '+str(len(db_meta)))
    return db_meta

""" This method serves as the entry point for the meta database. """
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

""" This is the main and starting point of the script. """
if __name__ == '__main__':
    release = get_current_release()
    if release:
        version = release.split('\t')
        print('Rapid release version = '+str(version[2])+ ' and release date is '+str(version[3]))
        rapid_meta = get_rapid_dbs(version[2])
        existing_assemblies = get_pending_genebuilds_for_update()
        update_genebuild_status(existing_assemblies,rapid_meta,version[2],version[3],os.getenv('GBS1'),os.getenv('REG_DB'),os.getenv('GBUSER'),int(os.getenv('GBP1')),os.getenv('GBPASS'))

    else:
        print('Could not determine current release version. Check that access to production meta server is good.\n')
