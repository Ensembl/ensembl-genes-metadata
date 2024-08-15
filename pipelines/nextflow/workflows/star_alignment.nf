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

if (!params.bam2cram) {
    exit 1, "Undefined --params.transcriptomic_dbname parameter. Please provide the server host for the db connection"
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
    log.info ' nextflow -C ensembl-genes-metadata/nextflow_star.config run nextflow/workflows/star_alignment.nf -entry STAR_ALIGNMENT  '
    log.info ''
    log.info 'Options:'
    log.info '  --bam2cram STR                   Oprion to convert BAM file to CRAM format  '
    log.info '  --outDir STR                 Output directory. Default is workDir'
    log.info '  --csvFile STR                Path for the csv containing the db name' 
    exit 1
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { FETCH_GENOME } from '../modules/star_alignment/fetch_genome.nf'
include { INDEX_GENOME } from '../modules/star_alignment/index_genome.nf'
include { DOWNLOAD_PAIRED_FASTQS } from '../modules/star_alignment/download_paired_fastqs.nf'
include { RUN_STAR } from '../modules/star_alignment/run_star.nf'
include { INDEX_BAM } from '../modules/star_alignment/index_bam.nf'
include { CONVERT_BAM_TO_CRAM } from '../modules/star_alignment/convert_bam_to_cram.nf'
include { INDEX_CRAM } from '../modules/star_alignment/index_cram.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow STAR_ALIGNMENT {
    def data = Channel.fromPath(params.csvFile, type: 'file', checkIfExists: true)
                .splitCsv(sep:',', header:true)
                .map { row -> [taxon_id:row.get('taxon_id'), gca:row.get('gca'), \
                run_accession:row.get('run_accession'), pair1:row.get('pair1'),
                pair2:row.get('pair2'),md5_1:row.get('md5_1'),md5_2:row.get('md5_2')]}
    data.each { dataRow -> dataRow.view() }    

    def genomeAndDataToAlign = FETCH_GENOME(data.flatten())
    def genomeIndexAndDataToAlign = INDEX_GENOME(genomeAndDataToAlign)
    def pairedFastqFiles=DOWNLOAD_PAIRED_FASTQS(genomeIndexAndDataToAlign)
    def starOutput = RUN_STAR(pairedFastqFiles)  
    def indexBamOutput = INDEX_BAM (starOutput)   
    if (params.bam2cram){
        def cramFile = CONVERT_BAM_TO_CRAM (indexBamOutput)
        def indexCramFile =  INDEX_CRAM (cramFile)
    }   

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
}

def deleteRecursively(Path path) {
    if (Files.isDirectory(path)) {
        Files.newDirectoryStream(path).each { subPath ->
            deleteRecursively(subPath)
        }
    }
    Files.delete(path)
}
