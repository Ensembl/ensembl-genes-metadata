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
include { checkRunStatus } from '../modules/utils.nf'
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
    //def run_Id = []
    //def bb=pairedFastqFilesMetadata
    //bb.each{ d-> d.view()}
    //d -> "GCA: ${d.gca}, Taxon ID: ${d.taxon_id}, run: ${d.run_accession},${pair1},${pair2},${dataFileQuery}"}
    def fastqcOutput = RUN_FASTQC(pairedFastqFilesMetadata)
    def processedFastQCOutput = fastqcOutput.map { result ->
//    getRunAccession.view { result ->
        def (taxon_id, gca, run_accession, pair1, pair2, dataFileQuery, fastqcPath) = result
        println "results: ${result}"
        println "Run Accession: ${run_accession}"
        // You can use run_accession here to get runId
        def run_Id = getDataFromTable("run_id", "run", "run_accession", run_accession)[0].run_id
        //def run_Id = getRunTable(run_accession, 'run_id')
        log.info "Run ID: ${run_Id}"
        return tuple(taxon_id, gca, run_accession, pair1, pair2, dataFileQuery, fastqcPath, run_Id.toString())
    }
    
    //def runId = getRunId(runAccession)
    def (runAccessionData, pairedFastqFiles, insertIntoDataFile) = PROCESS_FASTQC_OUTPUT(processedFastQCOutput)
    def updateValue = "False"
    def (runAccessionData_output,insertIntoDataFileQuery) = BUILD_QUERY(runAccessionData, insertIntoDataFile, updateValue)
    insertIntoDataFileQuery.subscribe { line ->
        def queriesArray = line.toString().split(";")
        setMetaDataRecord(queriesArray[0]+';')
        setMetaDataRecord(queriesArray[1]+';')
    } 

    def runAccessionData_QCstatus = runAccessionData_output.map { result ->
        def (taxon_id, gca, run_accession) = result
        def run_Id = getDataFromTable("run_id", "run", "run_accession", run_accession)[0].run_id
        checkRunStatus(runId)
        return tuple(taxon_id, gca, run_accession)
    }

    //if (qc_status == 'QC_PASS') {
    def subsamplingOutput = SUBSAMPLE_FASTQ_FILES(runAccessionData_QCstatus, pairedFastqFiles)

    emit: subsamplingOutputMetadata : subsamplingOutput
    /* NOT NEEDED FORE NOW
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

