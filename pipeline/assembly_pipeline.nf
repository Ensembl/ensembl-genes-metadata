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
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SET_DATE } from '../modules/set_date.nf'
include { FETCH_GCA } from '../modules/fetch_gca.nf'
include { PARSE_METADATA } from '../modules/parse_metadata.nf'
include { WRITE2DB_ASSEMBLY } from '../modules/write2db_assembly.nf'
include { UPDATE_KEYS_METADATA } from '../modules/update_keys_metadata.nf'
include { WRITE2DB_METADATA } from '../modules/write2db_metadata.nf'
include { SPECIES_CHECKER } from '../modules/species_checker.nf'
include { WRITE2DB_SPECIES } from '../modules/write2db_species.nf'
include { GET_TOLID } from '../modules/get_tolid.nf'
include { WRITE2DB_TOLID } from '../modules/write2db_tolid.nf'
include { CUSTOM_GROUP } from '../modules/custom_groups.nf'
include { WRITE2DB_GROUP } from '../modules/write2db_group.nf'
include { REPORT } from '../modules/report.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
WORKFLOW: REGISTER NEW ASSEMBLIES IN DB
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    // help
    if (params.help) {
    log.info ''
    log.info 'Pipeline to run Assembly metadata pipeline'
    log.info '-------------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'nextflow -C ensembl-genes-metadata/conf/assembly_pipeline.config \
                run ensembl-genes-metadata/pipeline/assembly_pipeline.nf \
                --enscode $ENSCODE --output_dir <OutDir>--taxon <taxon>'
    log.info ''
    log.info 'Options:'
    log.info '  --enscode STR               ENSCODE directory path'
    log.info '  --output_dir STR            Output directory path'
    log.info '  --date STR                  Custom date to retrieve assemblies (optional)'
    log.info '  --full_screen BOOLEAN       Run full screen mode, it will retrieve assemblies since 2019'
    log.info '  --taxon INT                 NCBI taxon id. Default is 2759'
    log.info '  --add_gca BOOLEAN           If option add_gca is set to true, it will use as input the GCA list provided by --gca_list'
    log.info '  --gca_list STR              GCA list file path. Requires --add_gca to be used as input'
    log.info '  --help BOOLEAN              Help option'
    }

    // print params
    params.each{ k, v -> println "params.${k.padRight(25)} = ${v}" }

    SET_DATE()
    def last_update = SET_DATE.out.splitText() {it -> it.trim()}

    FETCH_GCA(params.taxon, last_update)
    def gca = FETCH_GCA.out.splitText().map{it -> it.trim()}

    def parse_metadata_out = PARSE_METADATA(gca)

    def write2db_assembly_out = WRITE2DB_ASSEMBLY(parse_metadata_out)

    def update_keys_out = UPDATE_KEYS_METADATA(write2db_assembly_out)

    def write2db_metadata_out = WRITE2DB_METADATA(update_keys_out)

    def species_checker_out = SPECIES_CHECKER(write2db_metadata_out)

    def write2db_species_out = WRITE2DB_SPECIES(species_checker_out)

    def get_tolid_out = GET_TOLID(write2db_species_out)

    def write2db_tolid_out = WRITE2DB_TOLID(get_tolid_out)

    def custom_group_out = CUSTOM_GROUP(write2db_tolid_out)

    WRITE2DB_GROUP(custom_group_out)

    gca_list = WRITE2DB_GROUP.out.gca.map { it -> it.trim() }.collectFile(name: 'gca_list_to_report.txt', newLine: true).view()

    REPORT(gca_list, last_update)

}