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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES AND CONFIGURATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FETCH_ASSEMBLIES } from '../modules/fetch_assemblies.nf'
include { FETCH_METADATA } from '../modules/fetch_metadata.nf'
include { ASSEMBLY_STATUS } from '../modules/assembly_status.nf'
include { ASSEMBLY_REFSEQ } from '../modules/assembly_refseq.nf'
include { ASSEMBLY_METRICS } from '../modules/assembly_metrics.nf'
include { ASSEMBLY_NAME } from '../modules/assembly_name.nf'
include { BIOPROJECT } from '../modules/bioproject.nf'
include { TAXONOMY_CHECK } from '../modules/taxonomy_check.nf'
include { SPECIES_CHECKER } from '../modules/species_checker.nf'
include { WRITE2DB } from '../modules/write2db.nf'
include { TAXONOMY ; TAXONOMY as NEW_TAXONOMY } from '../modules/taxonomy.nf'
include { REPORT_UPDATE } from '../modules/report_update.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
WORKFLOW: REGISTER NEW ASSEMBLIES IN DB
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow ASSEMBLY_METADATA_UPDATE {
    // help
    if (params.help) {
    log.info"""
    ======================================================================
            Nextflow Pipeline to run Assembly metadata pipeline
    =======================================================================
    
    Usage: 
    nextflow -C ensembl-genes-metadata/conf/assembly_pipeline.config \
                run ensembl-genes-metadata/pipeline/assembly_pipeline.nf \
                --enscode $ENSCODE --output_dir <OutDir> --taxon <taxon>

    Required arguments:
    --output_dir STR            Output directory path

    Optional arguments:
    --screen_date STR           Custom date to retrieve assemblies and attempt update 
    --full_screen BOOLEAN       Run full screen mode, it will retrieve assemblies since 2019
    --gca_list STR              GCA list file path. Requires --add_gca to be used as input
    --help BOOLEAN              Help option
    """.stripIndent()
    }

    // print params
    params.each { k, v -> println "params.${k.padRight(25)} = ${v}" }

    main:

    FETCH_ASSEMBLIES(params.screen_date)
    def gca = FETCH_ASSEMBLIES.out.splitText().map { it -> it.trim() }

    def fetch_metadata_out = FETCH_METADATA(gca)

    ASSEMBLY_STATUS(fetch_metadata_out)
    ASSEMBLY_REFSEQ(fetch_metadata_out)
    ASSEMBLY_METRICS(fetch_metadata_out)
    ASSEMBLY_NAME(fetch_metadata_out)
    BIOPROJECT(fetch_metadata_out)

    def taxonomy_check_out = TAXONOMY_CHECK(fetch_metadata_out)

    taxonomy_check_out
        .branch { tuple ->
            def (gca_value, attempt_update, metadata_json, old_taxon_id, new_taxon_id, status) = tuple
            pass: status == 'pass'
                return [gca_value, attempt_update, metadata_json]
            failed: status == 'fail'
                return [gca_value, attempt_update, metadata_json, old_taxon_id, new_taxon_id]
        }
        .set { taxonomy_check_results }


    TAXONOMY(taxonomy_check_results.pass)

    def species_checker_out = SPECIES_CHECKER(taxonomy_check_results.failed)
    WRITE2DB(species_checker_out)
    NEW_TAXONOMY(WRITE2DB.out.to_taxonomy)

    def all_output = ASSEMBLY_STATUS.out.mix(ASSEMBLY_REFSEQ.out, ASSEMBLY_METRICS.out, ASSEMBLY_NAME.out, BIOPROJECT.out, TAXONOMY.out, NEW_TAXONOMY.out)
    .splitCsv()
    .map { row -> tuple(row[0].trim(), row[1].trim(), row[2].trim(), row[3].trim(), row[4].trim()) }
    .multiMap { item ->
        report: item
        tracking: item
    }
    
    if (params.slack_report) {
    REPORT_UPDATE(all_output.report)
    }

    all_output.tracking
    .collectFile(
    name: "${params.output_dir}/report_track.csv",
    seed: 'assembly,check_type,reporting,previous_value,new_value\n' ) { row -> row.join(',') + '\n' }

}

workflow.onComplete {
    log.info "Pipeline completed at: ${new Date().format('dd-MM-yyyy HH:mm:ss')}"
}