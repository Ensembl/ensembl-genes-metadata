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

import java.nio.file.*
import java.nio.file.attribute.BasicFileAttributes
import java.time.LocalDateTime
import java.time.format.DateTimeFormatter

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

//if (!params.enscode) {
//    exit 1, "Undefined --enscode parameter. Please provide the enscode path"
//}

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
    log.info ' nextflow -C ensembl-genes-metadata/nextflow.config run nextflow/workflows/short_read.nf -entry SHORT_READ  '
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
    exit 1
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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
    //taxon id present or not? if yes get all new short read data after this date if not add it for the first time
    def taxonomyResults= PROCESS_TAXONOMY_INFO(data)
    def fastqFilesMetadata  = PROCESS_RUN_ACCESSION_METADATA(taxonomyResults).collect()
    dd=fastqFilesMetadata
    dd.each{ d-> d.view()}
    //def fastQCMetadata = FASTQC_PROCESSING(fastqFilesMetadata).collect()
    //RUN_ALIGNMENT(fastQCMetadata)
    RUN_ALIGNMENT(fastqFilesMetadata)
    //if (params.cleanCache) {
        // Clean cache directories
    //    exec "rm -rf ${params.cacheDir}/*"
    //}
}  

workflow.onComplete {
    log.info "Pipeline completed at: ${new Date().format('dd-MM-yyyy HH:mm:ss')}"
    if (params.cleanOutputDir) {    
    try {
        def outDir = Paths.get(params.outDir)

        Files.newDirectoryStream(outDir, "*").each { path ->
        if (!path.toString().endsWith(".gz")) {
            deleteRecursively(path)
            }
        }
        log.info "Cleaning process completed successfully."
    } catch (Exception e) {
        log.error "Exception occurred while executing cleaning command: ${e.message}", e
    }
    }
    if (params.backupDB) {
        //def backupFilePath = "${params.outDir}/${params.transcriptomic_dbname}_backup.sql"
        // Generate a timestamp
        def timestamp = LocalDateTime.now().format(DateTimeFormatter.ofPattern("yyyyMMdd_HHmmss"))
        def backupFilePath = "${params.outDir}/${params.transcriptomic_dbname}_backup_${timestamp}.sql"
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
}

def deleteRecursively(Path path) {
    if (Files.isDirectory(path)) {
        Files.newDirectoryStream(path).each { subPath ->
            deleteRecursively(subPath)
        }
    }
    Files.delete(path)
    //println "Deleted: ${path}"
}
