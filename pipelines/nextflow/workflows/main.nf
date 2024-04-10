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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELP
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.help) {
    log.info"""
    =========================================
    Transcriptomic pipeline for short read evaluation and filtering
    =========================================
    Usage
    -----   

    The typical command for running the pipeline is as follows:

    .. code-block:: nextflow
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
params.outDir = params.outDir ?: config.outDir
params.assembly_accession = params.assembly_accession ?: config.assembly_accession

if (!params.taxon_id) {
    error "You must provide a value for taxon_id either via command line or config."
}

if (!params.assembly_accession) {
    error "You must provide a value for assembly accession either via command line or config."
}

if (!params.transcriptomic_dbhost) {
    exit 1, "Undefined --host parameter. Please provide the server host for the db connection"
}

if (!params.transcriptomic_dbport) {
    exit 1, "Undefined --port parameter. Please provide the server port for the db connection"
}
if (!params.transcriptomic_dbuser) {
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
    log.info '  --host STR                   Db host server '
    log.info '  --port INT                   Db port  '
    log.info '  --user STR                   Db user  '
    log.info '  --enscode STR                Enscode path '
    log.info '  --outDir STR                 Output directory. Default is workDir'
    log.info '  --taxon id INT               Taxon Id'
    log.info '  --assembly_accession STR     Assembly accession'
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
include { FETCH_GENOME } from '../modules/fetch_genome.nf'
//
// SUBWORKFLOW
//

include { PROCESS_TAXONOMY_INFO } from '../subworkflows/process_taxonomy_info.nf'
include { PROCESS_RUN_ACCESSION_METADATA } from '../subworkflows/process_run_accession_metadata.nf'
include { FASTQ_PROCESSING } from '../subworkflows/fastq_processing.nf'
include { RUN_ALIGNMENT } from '../subworkflows/run_alignment.nf'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
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
    genome_file=FETCH_GENOME(params.assembly_accession)
    //run_accession_list = key from PROCESS_TAXONOMY_INFO.runAccessionsPath.splitJson().flatten()
    run_accession_list = Channel.fromPath(PROCESS_TAXONOMY_INFO.runAccessionsPath).splitCsv()
    
    run_accession_list.subscribe { accession ->
        // Generate a job for each accession
        def processMetadata = PROCESS_RUN_ACCESSION_METADATA(params.taxon_id, accession)
        def processFastq = FASTQ_PROCESSING(processMetadata.taxon_id, accession, processMetadata.pairedFastqFiles)
        RUN_ALIGNMENT(params.assembly_accession,genome_file,processFastq.runAccessionFastqs)
    }
}  

workflow.onComplete {
    log.info "Pipeline completed at: ${new Date().format('dd-MM-yyyy HH:mm:ss')}"
}