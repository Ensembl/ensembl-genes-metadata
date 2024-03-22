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
========================================================================================
                         SteadyFlow - Steady Sate Transcription PIPELINE
========================================================================================
Steady Sate Analysis Pipeline. Started 2018-06-21.
 #### Homepage / Documentation
 https://github.com/Dowell-Lab/RNA-seq-Flow
 #### Authors
 Margaret Gruca <magr0763@colorado.edu>
========================================================================================
========================================================================================

Pipeline steps:

    1. Pre-processing sra/fastq
        # Proces taxonomy info
        # Download fastq/files (paired only)

    2. FastQC mapping and quality control
        # FastQC -- perform FastQC on fastq files and evaluate the status

    3. STAR alignment and mapping coverage evaluation
        # STAR - calculate mapping coverage for a read pair for a provided genome file 

*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELP
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.help) {
    log.info"""
    =========================================
    Transcriptomic pipeline for short read evaluation and filtering
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run main.nf -profile slurm --fastqs '/project/*_{R1,R2}*.fastq' --outdir '/project/'

    Required arguments:
        -profile                      Configuration profile to use. <base, slurm>
        --workdir                     Nextflow working directory where all intermediate files are saved.
        --taxon_id                    Input taxon id 
        --run_accession    specified as ['SRA..'] for now

    FastQC Options:

    STAR Options:

    
    Transcriptomic Db Connection:
        --dbname
        --host
        --port
        --user
        ---password

    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

includeConfig './nextflow.config'


params.taxon_id = params.taxon_id ?: config.taxon_id
params.genome_file = config.genome_file //it can be downloaded using assembly name and GCA
params.run_accession = = params.run_accession ?: config.run_accession  //TEMPORARY
params.outDir = params.outDir ?: config.outDir
params.assembly_name = params.assembly_name ?: config.assembly_name
params.assembly_accession = params.assembly_name ?: config.assembly_name

if (!params.taxon_id) {
    error "You must provide a value for taxon_id either via command line or config."
}

if (!params.gca) {
    error "You must provide a value for assembly accession either via command line or config."
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


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
    https://github.com/Ensembl/ensembl-genomio/blob/main/pipelines/nextflow/workflows/dumper_pipeline/main.nf
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/
//
include {                  } from '../modules/local/'

//
// SUBWORKFLOW
//

include { PROCESS_TAXONOMY_INFO } from '../subworkflows/process_taxonomy_info/process_taxonomy_info.nf'
include { PROCESS_RUN_ACCESSION_METADATA } from '../subworkflows/process_run_accession_metadata/process_run_accession_metadata.nf'
include { FASTQ_PROCESSING } from '../subworkflows/fastq_processing/fastq_processing.nf'
include { RUN_ALIGNMENT } from '../subworkflows/run_alignment/run_alignment.nf'
include { FETCH_GENOME } from '../modules/fetch_genome.nf'
//include { DUMP_SQL } from '../../subworkflows/dump_sql/main.nf'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW

    input taxon id, gca, assembly name, species name if needed for the genome file 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SHORT_READ {
    //jolly channel
    run_accession_list = Channel.empty()

    //taxon id present or not? if yes get all new short read data after this date if not add it for the first time
    //given taxon id, get list of run accession an 
    PROCESS_TAXONOMY_INFO (
        params.taxon_id,
    )
    genome_file=FETCH_GENOME(gca)
    //run_accession_list = key from PROCESS_TAXONOMY_INFO.runAccessionsPath.splitJson().flatten()
    run_accession_list = Channel.fromPath(PROCESS_TAXONOMY_INFO.runAccessionsPath).splitCsv()
    
    run_accession_list.subscribe { accession ->
        // Generate a job for each accession
        def processMetadata = PROCESS_RUN_ACCESSION_METADATA(params.taxon_id, accession)
        def processFastq = FASTQ_PROCESSING(processMetadata.taxon_id, accession, processMetadata.pairedFastqFiles)
        RUN_ALIGNMENT(params.gca,genome_file,processFastq.runAccessionFastqs)
    }
}  

workflow.onComplete {
    log.info "Pipeline completed at: ${new Date().format('dd-MM-yyyy HH:mm:ss')}"
}