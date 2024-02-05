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

if (!params.taxon_id) {
    error "You must provide a value for taxon_id either via command line or config."
}

if (!params.run_accession) {
    error "You must provide a value for taxon_id either via command line or config."
}


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

include { PROCESS_TAXONOMY_INFO } from '../../subworkflows/process_taxonomy_info/process_taxonomy_info.nf'
include { PROCESS_RUN_ACCESSION_METADATA } from '../../subworkflows/process_run_accession_metadata/process_run_accession_metadata.nf'
include { FASTQ_PROCESSING } from '../../subworkflows/fastq_processing/fastq_processing.nf'
include { RUN_ALIGNMENT } from '../../subworkflows/run_alignment/run_alignment.nf'
//include { DUMP_SQL } from '../../subworkflows/dump_sql/main.nf'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SHORT_READ {
    //jolly channel
    run_accession_list = Channel.empty()
    if (!params.run_accession ) {
    //given taxon id, gete list of run accession an 
    PROCESS_TAXONOMY_INFO (
        params.taxon_id,
        params.transcriptomic_dbname ,                
        params.transcriptomic_host,                  
        params.transcriptomic_port,   
        params.transcriptomic_user,    
        params.transcriptomic_password,

    )
    run_accession_list = run_accession_list(PROCESS_TAXONOMY_INFO.out.filtered_run_accessions).flatten()

    }
    else {
      run_accession_list = params.run_accession
    }
    //get and store metadat and download fastq paired files 
    paired_fastq= PROCESS_RUN_ACCESSION_METADATA (
      run_accession_list

    )

    csvData = Channel.fromPath(params.csvFile).splitCsv()

    // Get db name and its metadata
    db = csvData.flatten()
    db_meta = SPECIES_METADATA(db, params.outDir, params.project)
      .splitCsv(header: true)

    // Get the closest Busco dataset from the taxonomy classification stored in db meta table 
    db_dataset = BUSCO_DATASET(db_meta)
    
    // Run Busco in genome mode
    if (busco_mode.contains('genome')) {
        genome_data = FETCH_GENOME(db_dataset, params.cacheDir)
        busco_genome_output = BUSCO_GENOME_LINEAGE(genome_data)
        BUSCO_GENOME_OUTPUT(busco_genome_output, "genome", params.project)
        if (params.project == 'ensembl') {
          FASTA_GENOME_OUTPUT(genome_data, params.project, 'genome')
        }
    }
    
    // Run Busco in protein mode
    if (busco_mode.contains('protein')) {
        if (params.project == 'brc') {
            db_dataset = db_dataset.filter{ it[0].has_genes == "1" }
        }
        protein_data = FETCH_PROTEINS (db_dataset, params.cacheDir)
        busco_protein_output = BUSCO_PROTEIN_LINEAGE(protein_data)
        BUSCO_PROTEIN_OUTPUT(busco_protein_output, "protein", params.project)
        if (params.project == 'ensembl') {
          FASTA_PROTEIN_OUTPUT(protein_data, params.project, 'fasta')
        }
    }
}