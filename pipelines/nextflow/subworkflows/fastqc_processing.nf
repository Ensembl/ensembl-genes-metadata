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


includeConfig './pipelines/workflows/nextflow.config'
include { checkFastqc, checkOverrepresentedSequences  } from './pipelines/modules/utils.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { PROCESS_FASTQC_OUTPUT } from '../modules/fastqc_processing/process_fastqc_output.nf'
include { RUN_FASTQC } from '../modules/fastqc_processing/run_fastqc.nf'
include { BUILD_QUERY } from '../modules/store_metadata.nf'
include { SUBSAMPLE_FASTQ_FILES } from '../modules/fastqc_processing/subsample_fastq_files.nf'
include { ADAPTER_TRIMMING } from '../modules/fastqc_processing/adapter_trimming.nf'
include { CLEANING } from '../modules/cleaning.nf'

include { setMetaDataRecord } from '../modules/utils.nf'
include { getDataFromTable } from '../modules/utils.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FASTQC_PROCESSING{
    take:        
    val pairedFastqFilesMetadata 
    
    main:
    def fastqcOutput = RUN_FASTQC(pairedFastqFilesMetadata.flatten())
    def runAccession = fastqcOutput.run_accession
    def runId = getRunId(runAccession)
    def (runAccessionData, pairedFastqFiles, insertIntoDataFile) = PROCESS_FASTQC_OUTPUT(fastqcOutput, runId)
    def updateValue = "False"
    def (runAccessionData_1,insertIntoDataFileQuery) = BUILD_QUERY(runAccessionData, insertIntoDataFile, updateValue)
    insertIntoDataFileQuery.subscribe { line ->
        setMetaDataRecord(line.toString())
    }

    def qc_status = getDataFromTable(runAccessionData_1.run_accession, "run", "qc_status")

    emit:
    qc_status
    if (qc_status == 'QC_PASS'){
        subsampling = SUBSAMPLE_FASTQ_FILES(runAccessionData, pairedFastqFiles)
        def  overrepresented_sequences = checkOverrepresentedSequences(params.jdbcUrl, params.transcriptomic_dbuser, params.transcriptomic_dbpassword, run_accession)

        if (overrepresented_sequences == true){
          trimming = ADAPTER_TRIMMING(subsampling.subsampledFastqs.taxon_id, subsampling.run_accession, subsampling.subsampledFastqs)
          emit runAccessionFastqs(trimming.runTrimmedFastqs.taxon_id, trimming.runTrimmedFastqs.run_accession, trimming.runTrimmedFastqs.trimmedFastqFiles)
        } else {
          emit runAccessionFastqs(subsampling.sampledFastqFiles.taxon_id, subsampling.sampledFastqFiles.run_accession, subsampling.sampledFastqFiles.pairedFastqFiles)
        }
    } else if (qc_status == 'QC_FAIL'){
      CLEANING(taxon_id, run_accession)
    }
}
/*
    insertIntoRunQueryOutputs= insertIntoRunQuery.toString().split('\n')
    if (insertIntoRunQueryOutputs.size() == 1) {
        // Only one output, pass it directly
        setMetaDataRecord(insertIntoRunQuery.toString())
    } else {
        // Multiple outputs, loop through and pass each one
        insertIntoRunQueryOutputs.each { singleOutput ->
            setMetaDataRecord(singleOutput)
        }
        }
        */