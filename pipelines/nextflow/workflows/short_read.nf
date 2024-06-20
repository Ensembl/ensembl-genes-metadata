#!/usr/bin/env nextflow
/*
See the NOTICE file distributed with this work for additional information
regarding copyright ownership.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

nextflow.enable.dsl=2

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

if (!params.transcriptomic_dbname) {
    exit 1, "Undefined --params.transcriptomic_dbname parameter. Please provide the server host for the db connection"
}

if (!params.transcriptomic_dbhost) {
    exit 1, "Undefined --transcriptomic_dbhost parameter. Please provide the server host for the db connection"
}

if (!params.transcriptomic_dbport) {
    exit 1, "Undefined --transcriptomic_dbport parameter. Please provide the server port for the db connection"
}
if (!params.transcriptomic_dbuser) {
    exit 1, "Undefined --transcriptomic_dbuser parameter. Please provide the server user for the db connection"
}

if (!params.enscode) {
    exit 1, "Undefined --enscode parameter. Please provide the enscode path"
}

if (!params.cacheDir) {
    exit 1, "Undefined --cacheDir parameter. Please provide the cache dir directory's path"
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELP
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.help) {
    log.info ''
    log.info 'Pipeline to process short read data fo a given taxon id'
    log.info '-------------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '  nextflow -C ...'
    log.info ''
    log.info 'Options:'
    log.info '  --transcriptomic_dbname STR                   Db name '
    log.info '  --transcriptomic_dbhost STR                   Db host server '
    log.info '  --transcriptomic_dbport INT                   Db port  '
    log.info '  --transcriptomic_dbuser STR                   Db user  '
    log.info '  --transcriptomic_dbpassword STR                   Db password  '
    log.info '  --user_r STR                 Db user read_only'
    log.info '  --enscode STR                Enscode path '
    log.info '  --outDir STR                 Output directory. Default is workDir'
    log.info '  --csvFile STR                Path for the csv containing the db name' 
    log.info '  --bioperl STR                BioPerl path (optional)'
    exit 1
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
    https://github.com/Ensembl/ensembl-genomio/blob/main/pipelines/nextflow/workflows/dumper_pipeline/main.nf
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE
//
//
// SUBWORKFLOW
//

include { PROCESS_TAXONOMY_INFO } from '../subworkflows/process_taxonomy_info.nf'
include { PROCESS_RUN_ACCESSION_METADATA } from '../subworkflows/process_run_accession_metadata.nf'
include { FASTQC_PROCESSING } from '../subworkflows/fastqc_processing.nf'
include { RUN_ALIGNMENT } from '../subworkflows/run_alignment.nf'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SHORT_READ {
    def data = Channel.fromPath(params.csvFile, type: 'file', checkIfExists: true)
                .splitCsv(sep:',', header:true)
                .map { row -> [taxon_id:row.get('taxon_id'), gca:row.get('gca')]}
    data.each { dataRow -> dataRow.view() }            
    //data1=data
    //data1.each{ d-> d.view()}
    //taxon id present or not? if yes get all new short read data after this date if not add it for the first time
    //given taxon id, get list of run accession an 
    def taxonomyResults= PROCESS_TAXONOMY_INFO(data)
    //rr=runAccessionList
    //rr.each { d -> "Taxon ID: ${d.taxon_id}, GCA: ${d.gca}, run accession: ${d.run_accession}"} 
    def fastqFilesMetadata  = PROCESS_RUN_ACCESSION_METADATA(taxonomyResults).collect()
    dd=fastqFilesMetadata
    dd.each{ d-> d.view()}
    //def fastqMetadata = []
    //fastqMetadata.add(fastqFilesMetadata)
    //FASTQC_PROCESSING(fastqFilesMetadata)
    def fastQCMetadata = FASTQC_PROCESSING(fastqFilesMetadata).collect()
    
    RUN_ALIGNMENT(fastQCMetadata)
    //Channel.of(runAccessionList).groupTuple(taxon_id).collect()
    //runAccessionList.view()
//    if (params.backupDB) {
 //       exec """pg_dump -h ${params.transcriptomic_dbhost} -p ${params.transcriptomic_dbport} -U ${params.transcriptomic_dbuser} ${params.transcriptomic_dbname} | gzip > ${params.outDir}/${params.transcriptomic_dbname}_backup.sql.gz"""
   //     }
   /* 
    run_accession_list.subscribe { accession ->
        // Generate a job for each accession
        def processMetadata = PROCESS_RUN_ACCESSION_METADATA(params.taxon_id, accession)
        def processFastq = FASTQ_PROCESSING(processMetadata.taxon_id, accession, processMetadata.pairedFastqFiles)
        RUN_ALIGNMENT(params.assembly_accession,genome_file,processFastq.runAccessionFastqs)
    }
*/
    if (params.cleanCache) {
        // Clean cache directories
        exec "rm -rf ${params.cacheDir}/*"
    }
}  

workflow.onComplete {
    log.info "Pipeline completed at: ${new Date().format('dd-MM-yyyy HH:mm:ss')}"
    def cleaningCommand= ['rm','-rf', params.outDir+'/*']

    log.info "Executing cleaning command: ${cleaningCommand.join(' ')}"

    def cleaningProcess = new ProcessBuilder(cleaningCommand).start()
    cleaningProcess.waitFor()
    if (params.backupDB) {
        def backupFilePath = "${params.outDir}/${params.transcriptomic_dbname}_backup.sql"
def gzipFilePath = "${backupFilePath}.gz"

// Define the database backup command as a list of arguments
def mysqldumpCommand = [
    'mysqldump',
    '-h', params.transcriptomic_dbhost,
    '-P', params.transcriptomic_dbport.toString(),
    '-u', params.transcriptomic_dbuser,
    '-p' + params.transcriptomic_dbpassword,
    params.transcriptomic_dbname
]

log.info "Executing database backup command: ${mysqldumpCommand.join(' ')}"

def processBuilder = new ProcessBuilder(mysqldumpCommand)
processBuilder.redirectOutput(new File(backupFilePath))
processBuilder.redirectErrorStream(true) // Combine stdout and stderr

def mysqldumpProcess = processBuilder.start()
mysqldumpProcess.waitFor()

if (mysqldumpProcess.exitValue() != 0) {
    log.error "Database backup failed. See error output for details."
    mysqldumpProcess.inputStream.eachLine { line -> log.error line }
} else {
    log.info "Database backup completed successfully. Proceeding to gzip the backup file."
    
    def gzipCommand = [
        'gzip', '-f', backupFilePath.toString()
    ]

    log.info "Executing gzip command: ${gzipCommand.join(' ')}"
    
    def gzipProcess = new ProcessBuilder(gzipCommand).start()
    gzipProcess.waitFor()
    
    if (gzipProcess.exitValue() != 0) {
        log.error "Gzip failed. See error output for details."
        gzipProcess.inputStream.eachLine { line -> log.error line }
    } else {
        log.info "Gzip completed successfully. Backup file: ${gzipFilePath}"
    }
}

            }
            //mysqldump -h mysql-ens-genebuild-prod-1 -P 4527 -u ensadmin -pensembl transcriptomic_fra | gzip > transcriptomic_fra_backup.sql.gz
              
}
