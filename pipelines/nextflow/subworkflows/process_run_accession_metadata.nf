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

/*
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   VALIDATE INPUTS
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (!params.host) {
  exit 1, "Undefined --host parameter. Please provide the server host for the db connection"
}

if (!params.port) {
  exit 1, "Undefined --port parameter. Please provide the server port for the db connection"
}
if (!params.user) {
  exit 1, "Undefined --user parameter. Please provide the server user for the db connection"
}

if (!params.enscode) {
  exit 1, "Undefined --enscode parameter. Please provide the enscode path"
}
if (!params.outDir) {
  exit 1, "Undefined --outDir parameter. Please provide the output directory's path"
}
if (!params.cacheDir) {
  exit 1, "Undefined --cacheDir parameter. Please provide the cache dir directory's path"
}
if (!params.mode) {
  exit 1, "Undefined --mode parameter. Please define Busco running mode"
}

if (params.csvFile) {
    csvFile = file(params.csvFile, checkIfExists: true)
} else {
    exit 1, 'CSV file not specified!'
}

busco_mode = []
if (params.mode instanceof java.lang.String) {
  busco_mode = [params.mode]
}
else {
  busco_mode = params.mode
}

acceptable_projects = ['ensembl', 'brc']
if (!acceptable_projects.contains(params.project)) {
    exit 1, 'Invalid project name'
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELP
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.help) {
  log.info ''
  log.info 'Pipeline to run Busco score in protein and/or genome mode'
  log.info '-------------------------------------------------------'
  log.info ''
  log.info 'Usage: '
  log.info '  nextflow -C ensembl-genes-nf/nextflow.config run ensembl-genes-nf/iworkflows/busco_pipeline.nf --enscode --csvFile --mode'
  log.info ''
  log.info 'Options:'
  log.info '  --host STR                   Db host server '
  log.info '  --port INT                   Db port  '
  log.info '  --user STR                   Db user  '
  log.info '  --enscode STR                Enscode path '
  log.info '  --outDir STR                 Output directory. Default is workDir'
  log.info '  --csvFile STR                Path for the csv containing the db name'
  log.info '  --mode STR                   Busco mode: genome or protein, default is to run both'
  log.info '  --bioperl STR                BioPerl path (optional)'
  log.info '  --project STR                Project, for the formatting of the output ("ensembl" or "brc")'
  exit 1
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { STORE_METADATA } from '../modules/process_run_accession_metadata/store_metadata.nf'
include { DOWNLOAD_PAIRED_FASTQS } from '../modules/process_run_accession_metadata/download_paired_fastqs.nf'
include { GET_RUN_ACCESSION_METADATA } from '../modules/process_run_accession_metadata/get_run_accession_metadata.nf'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PROCESS_RUN_ACCESSION_METADATA {
    take:
    val taxon_id
    val run_accession          

    main:

    paired_fastq_files_path = Channel.empty()
    //it is an insert but we need to split the value 
    //so it might be  first function that split the values and another function INSERT
    //emit fasta_paired_files
    run_accession_metadata = GET_RUN_ACCESSION_METADATA(run_accession)
    STORE_METADATA(run_accession_metadata.runAccessionsMetadataPath)
    paired_fastq_files_path=DOWNLOAD_PAIRED_FASTQ(taxon_id, run_accession)

    emit:
    pairedFastqFiles            = paired_fastq_files_path                  // channel: [path1, path2]
}