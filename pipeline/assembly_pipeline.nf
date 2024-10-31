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


if (!params.output_dir) {
    exit 1, "Undefined --output_dir parameter. Please provide the output directory's path"
}

if (!params.enscode) {
    exit 1, "Undefined --enscode parameter. Please provide the enscode path"
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELP
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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
    log.info '  --help BOOLEAN              Help option'
}

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

// print params
params.each{ k, v -> println "params.${k.padRight(25)} = ${v}" }

workflow {

    SET_DATE()
    def last_update = SET_DATE.out.splitText() {it -> it.trim()}
    
    FETCH_GCA(params.taxon, last_update)

    gca = FETCH_GCA.out.splitText().map{it -> it.trim()}

    PARSE_METADATA(gca)
    
    WRITE2DB_ASSEMBLY(PARSE_METADATA.out.gca, PARSE_METADATA.out.assembly, PARSE_METADATA.out.metadata_tmp, PARSE_METADATA.out.species_tmp)
    
    UPDATE_KEYS_METADATA(WRITE2DB_ASSEMBLY.out.gca, WRITE2DB_ASSEMBLY.out.metadata_tmp, WRITE2DB_ASSEMBLY.out.last_id, WRITE2DB_ASSEMBLY.out.species_tmp)
    
    WRITE2DB_METADATA(UPDATE_KEYS_METADATA.out.gca, UPDATE_KEYS_METADATA.out.metadata, UPDATE_KEYS_METADATA.out.species_tmp)
    
    SPECIES_CHECKER(WRITE2DB_METADATA.out.gca, WRITE2DB_METADATA.out.species_tmp)
    
    WRITE2DB_SPECIES(SPECIES_CHECKER.out.gca, SPECIES_CHECKER.out.species)
    
    GET_TOLID(WRITE2DB_SPECIES.out.gca)
    
    WRITE2DB_TOLID(GET_TOLID.out.gca, GET_TOLID.out.tolid)
    
    CUSTOM_GROUP(WRITE2DB_TOLID.out.gca)
    
    WRITE2DB_GROUP(CUSTOM_GROUP.out.gca, CUSTOM_GROUP.out.group)

    gca_list = WRITE2DB_GROUP.out.gca.map { it -> it.trim() }.collectFile(name: 'gca_list_to_report.txt', newLine: true)

    REPORT(gca_list, last_update)

}