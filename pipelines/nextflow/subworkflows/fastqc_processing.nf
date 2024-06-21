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

import java.time.LocalDateTime
import java.time.format.DateTimeFormatter


//include { checkFastqc  } from './pipelines/modules/utils.nf'
//include { checkOverrepresentedSequences  } from './pipelines/modules/utils.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { PROCESS_FASTQC_OUTPUT } from '../modules/fastqc_processing/process_fastqc_output.nf'
include { RUN_FASTQC } from '../modules/fastqc_processing/run_fastqc.nf'
include { BUILD_QUERY } from '../modules/build_query.nf'
include { SUBSAMPLE_FASTQ_FILES } from '../modules/fastqc_processing/subsample_fastq_files.nf'
include { ADAPTER_TRIMMING } from '../modules/fastqc_processing/adapter_trimming.nf'
include { CLEANING } from '../modules/cleaning.nf'

include { setMetaDataRecord } from '../modules/utils.nf'
include { checkRunFastQCStatus } from '../modules/utils.nf'
include { checkOverrepresentedSequences } from '../modules/utils.nf'
include { getDataFromTable } from '../modules/utils.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: FAST QC PROCESSING AND SUBSAMPLING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/ 

workflow FASTQC_PROCESSING{
    take:        
    pairedFastqFilesMetadata 
    
    main:
    def bb=pairedFastqFilesMetadata
    bb.flatten().view {
    d -> "GCA: ${d.gca}, Taxon ID: ${d.taxon_id}, run: ${d.run_accession}, files: ${d.pair1},${d.pair2},${d.dataFileQuery}"}

    def fastqcOutput = RUN_FASTQC(pairedFastqFilesMetadata.flatten())
    def processedFastQCOutput = fastqcOutput.map { result ->
        def (taxon_id, gca, run_accession, dataFileQuery, fastqcPath) = result
        log.info("results: ${result}")
        log.info("Run Accession: ${run_accession}")

        def run_Id = getDataFromTable("run_id", "run", "run_accession", run_accession)[0].run_id
        log.info "Run ID: ${run_Id}"
        return tuple(taxon_id, gca, run_accession, dataFileQuery, fastqcPath, run_Id.toString())
    }
    
    def (runAccessionData, insertIntoDataFile) = PROCESS_FASTQC_OUTPUT(processedFastQCOutput)
    def updateValue = "False"
    def (runAccessionData_output,insertIntoDataFileQuery) = BUILD_QUERY(runAccessionData, insertIntoDataFile, updateValue)
    insertIntoDataFileQuery.subscribe { line ->
        log.info("insertIntoDataFileQuery.subscribe ${line}")
        def queriesArray = line.toString().split(";")
        queriesArray.eachWithIndex { query, index ->
            // Trim the query to remove any leading/trailing whitespace
            query = query.trim()
            // Check if the query is not empty
            if (query) {
                log.info("queriesArray[${index}] ${query};")
                setMetaDataRecord(query.toString() + ";")
            }
        }
    }

    def runAccessionData_QCstatus = runAccessionData_output.map { result ->
        def (taxon_id, gca, run_accession) = result
        log.info("runAccessionData_QCstatus ${run_accession}")
        def run_Id = getDataFromTable("run_id", "run", "run_accession", run_accession)[0].run_id.toString()
        checkRunFastQCStatus(run_Id)
        return tuple(taxon_id, gca, run_accession)
    }
    //if (qc_status == 'QC_PASS') {
    def subsamplingOutput = SUBSAMPLE_FASTQ_FILES(runAccessionData_QCstatus)
    def subsampling_Output = subsamplingOutput.map { subsampling ->
        def (taxon_id, gca, run_accession, sub1, sub2) = subsampling
        return [taxon_id:taxon_id, gca:gca, run_accession:run_accession, pair1:sub1, pair2:sub2]
    }    
    emit: subsamplingOutputMetadata = subsampling_Output
    
    /* NOT NEEDED FOR NOW
    def processedSamplingOutput = subsamplingOutput.map { result ->
        def (taxon_id, gca, run_accession, subPair1, subPair2) = result
        println "results: ${result}"
        println "Run Accession: ${run_accession}"
        
        def overrepresented_sequences = checkOverrepresentedSequences(run_accession)
        log.info "overrepresented_sequences: ${overrepresented_sequences}"
        
        return [taxon_id, gca, run_accession, subPair1, subPair2, overrepresented_sequences]
    }
    
    // Separate outputs based on the overrepresented_sequences status
    def subsampleWithOverrepresented = processedSamplingOutput.findAll { it[5] == true }
    def subsampleWithoutOverrepresented = processedSamplingOutput.findAll { it[5] == false }
    // overrepresented_sequences absent
    if (!subsampleWithoutOverrepresented.isEmpty()) {
        emit subsamplingOutputMetadata: subsampleWithoutOverrepresented
    }
    //overrepresented_sequences present, trimming needed
    if (!subsampleWithOverrepresented.isEmpty()) {
        def trimmingOutput = ADAPTER_TRIMMING(subsampleWithOverrepresented)
        emit trimmingOutputMetadata: trimmingOutput
    }
    
    //} else if (qc_status == 'QC_FAIL'){
    //  CLEANING(runAccessionData_output)
    //}
    */
    
}

