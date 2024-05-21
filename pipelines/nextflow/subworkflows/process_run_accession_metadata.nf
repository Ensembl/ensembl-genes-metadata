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


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { STORE_METADATA as STORE_RUN_METADATA } from '../modules/store_metadata.nf'
include { STORE_METADATA as STORE_STUDY_METADATA} from '../modules/store_metadata.nf'
include { DOWNLOAD_PAIRED_FASTQS } from '../modules/process_run_accession_metadata/download_paired_fastqs.nf'
include { GET_RUN_ACCESSION_METADATA } from '../modules/process_run_accession_metadata/get_run_accession_metadata.nf'

include { setMetaDataRecord } from '../modules/utils.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PROCESS_RUN_ACCESSION_METADATA {
    take:
    //taxon_id
    //run_accession          
    transcriptomic_meta
    //tuple val(taxon_id), val(gca), val(run_accession)
    main:
    def db_meta1=transcriptomic_meta
    db_meta1.flatten().view { d -> "GCA: ${d.gca}, Taxon ID: ${d.taxon_id}, run: ${d.run_accession}"}
    //it is an insert but we need to split the value 
    //so it might be  first function that split the values and another function INSERT
    //emit fasta_paired_files
    def(runAccessionMedatadata, insertIntoRun, insertIntoStudy, queryDataFile) = GET_RUN_ACCESSION_METADATA(transcriptomic_meta.flatten())
    def updateValue = "False"
    def (runAccessionMedatadata_1,insertIntoRunQuery) = STORE_RUN_METADATA(runAccessionMedatadata, insertIntoRun, updateValue)
    def data1=insertIntoRunQuery
    data1.flatten().view { d -> "query ${d}"}
    insertIntoRunQuery.subscribe { line ->

    setMetaDataRecord(line.toString())
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
    def (runAccessionMedatadata_2, insertIntoStudyQuery) = STORE_STUDY_METADATA(runAccessionMedatadata_1, insertIntoStudy, updateValue)
    setMetaDataRecord(insertIntoStudyQuery.toString())
    /*
    def insertIntoStudyQueryOutputs = insertIntoStudyQuery.toString().split('\n')
    if (insertIntoStudyQueryOutputs.size() == 1) {
        // Only one output, pass it directly
        setMetaDataRecord(insertIntoStudyQuery.toString())
    } else {
        // Multiple outputs, loop through and pass each one
        insertIntoStudyQueryOutputs.each { singleOutput ->
            setMetaDataRecord(singleOutput)
        }
        }
        */
    pairedFastqFiles=DOWNLOAD_PAIRED_FASTQS(runAccessionMedatadata_2)

    emit:
    insertIntoDataFile =   queryDataFile   // path
    pairedFastqFilesMetadata =   pairedFastqFiles      // channel: [taxon_id, gca, run_accession, path1, path2]
}
