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
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RUN_STAR } from '../modules/transcriptomic_alignment/run_star.nf'
include { EXTRACT_UNIQUELY_MAPPED_READS_PERCENTAGE } from '../modules/transcriptomic_alignment/extract_uniquely_mapped_reads_percentage.nf'
include { STORE_METADATA } from '../modules/store_metadata.nf'
include { CLEANING } from '../modules/cleaning.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    take:
    val gca
    val genome_file
    tuple val(taxon_id), val(run_accession), val(pairedFastqFiles)
    set fastqFile1, fastqFile2 from pairedFastqFiles

    main:

      genome_file=FETCH_GENOME(params.assembly_accession)
    star_output = RUN_STAR(genome_file.gca, genome_file.genome_file, fastqFile1, fastqFile2 )
    start_stats_json = EXTRACT_UNIQUELY_MAPPED_READS_PERCENTAGE(star_output.log_final_out, run_accession, gca)
    STORE_METADATA(start_stats_json.star_metadata)
    updateFastqcStatus(params.jdbcUrl, params.transcriptomic_dbuser, params.transcriptomic_dbpassword, run_accession)
    CLEANING(taxon_id, run_accession)

}

