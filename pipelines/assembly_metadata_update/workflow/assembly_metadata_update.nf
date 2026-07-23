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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES AND CONFIGURATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FETCH_ASSEMBLIES } from '../modules/fetch_assemblies.nf'
include { INTEGRITY_CHECKER } from '../modules/integrity_checker.nf'
include { INTEGRITY_TAXONOMY } from '../modules/integrity_taxonomy.nf'
include { INTEGRITY_WRITE2DB } from '../modules/integrity_write2db.nf'
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

    // screen_date is unused (and left unset) in --gca_input mode; substitute a
    // placeholder so a null value is never passed into the process input.
    def screen_date_value = params.gca_input ? 'NA' : params.screen_date
    FETCH_ASSEMBLIES(screen_date_value)
    def gca = FETCH_ASSEMBLIES.out.splitText().map { it -> it.trim() }

    INTEGRITY_CHECKER(gca)

    INTEGRITY_CHECKER.out
        .map { gca_value, stdout ->
            def line = stdout.trim()
            def parts = line.split(',')
            def status = parts[0].trim()
            def accession = parts[1].trim()
            return [gca_value, status, accession]
        }
        .branch { tuple ->
            def (gca_value, status, accession) = tuple
            correct: status == 'correct'
            return gca_value
            taxonomy_update: status == 'taxonomy_update'
            return [gca_value, accession]
            deleted: status == 'delete'
            return gca_value
            check: status == 'check'
            return gca_value
        }
        .set { integrity_check_results }

    INTEGRITY_TAXONOMY(integrity_check_results.taxonomy_update)
    INTEGRITY_WRITE2DB(INTEGRITY_TAXONOMY.out)

    def gca_accession = integrity_check_results.correct.mix(INTEGRITY_WRITE2DB.out.gca_to_update)

    def fetch_metadata_out = FETCH_METADATA(gca_accession)

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

    def all_output = ASSEMBLY_STATUS.out
        .mix(ASSEMBLY_REFSEQ.out, ASSEMBLY_METRICS.out, ASSEMBLY_NAME.out, BIOPROJECT.out, TAXONOMY.out, NEW_TAXONOMY.out)
        .splitCsv()
        .map { row -> tuple(row[0].trim(), row[1].trim(), row[2].trim(), row[3].trim(), row[4].trim()) }
        .multiMap { item ->
            report: item
            tracking: item
        }

    if (params.slack_report) {
        REPORT_UPDATE(all_output.report)
    }

    integrity_check_results.deleted.collectFile(
        name: "${params.output_dir}/deleted_GCAS_to_add.csv"
    ) { gca_value -> "${gca_value}\n" }

    integrity_check_results.check.collectFile(
        name: "${params.output_dir}/to_manually_check_GCAS.csv"
    ) { gca_value -> "${gca_value}\n" }


    all_output.tracking.collectFile(
        name: "${params.output_dir}/report_track.csv",
        seed: 'assembly,check_type,reporting,previous_value,new_value\n',
    ) { row -> row.join(',') + '\n' }

    workflow.onComplete {
        log.info("Pipeline completed at: ${new Date().format('dd-MM-yyyy HH:mm:ss')}")
    }
}
