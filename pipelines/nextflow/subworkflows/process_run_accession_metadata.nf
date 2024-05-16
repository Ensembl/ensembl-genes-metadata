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

    //it is an insert but we need to split the value 
    //so it might be  first function that split the values and another function INSERT
    //emit fasta_paired_files
    def(runAccessionMedatadata, insertIntoRun, insertIntoStudy, queryDataFile) = GET_RUN_ACCESSION_METADATA(transcriptomic_meta.flatten())
    STORE_RUN_METADATA(insertIntoRun)
    STORE_STUDY_METADATA(insertIntoStudy)
    pairedFastqFiles=DOWNLOAD_PAIRED_FASTQ(runAccessionMedatadata)

    emit:
    queryDataFile = insertIntoDataFile  // path
    pairedFastqFiles            = pairedFastqFilesMetadata              // channel: [taxon_id, gca, run_accession, path1, path2]
}
