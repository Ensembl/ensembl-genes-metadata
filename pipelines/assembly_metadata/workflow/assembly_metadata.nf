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
nextflow.enable.dsl = 2
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES AND CONFIGURATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SET_DATE } from '../modules/set_date.nf'
include { FETCH_GCA } from '../modules/fetch_gca.nf'
include { PARSE_METADATA } from '../modules/parse_metadata.nf'
include { UPDATE_KEYS_METADATA } from '../modules/update_keys_metadata.nf'
include { SPECIES_CHECKER } from '../modules/species_checker.nf'
include { GET_TOLID } from '../modules/get_tolid.nf'
include { REPORT } from '../modules/report.nf'
include {
    WRITE2DB as WRITE2DB_ASSEMBLY ;
    WRITE2DB as WRITE2DB_METADATA ;
    WRITE2DB as WRITE2DB_SPECIES ;
    WRITE2DB as WRITE2DB_TOLID
} from '../modules/write2db.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
WORKFLOW: REGISTER NEW ASSEMBLIES IN DB
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ASSEMBLY_METADATA {

    def write2db_script = file("${projectDir}/../../src/python/write2db.py")
    def species_checker_script = file("${projectDir}/../../src/python/species_checker.py")

    SET_DATE()
    def last_update = SET_DATE.out.splitText { it -> it.trim() }

    FETCH_GCA(params.taxon, last_update)
    def gca = FETCH_GCA.out.gca.splitText().map { it -> it.trim() }

    def parse_metadata_out = PARSE_METADATA(gca)

    def write2db_assembly_in = parse_metadata_out.map { gca_value, assembly, metadata_tmp, species_tmp -> tuple(gca_value, assembly, [metadata_tmp, species_tmp]) }
    WRITE2DB_ASSEMBLY(write2db_assembly_in, write2db_script, false)
    def write2db_assembly_out = WRITE2DB_ASSEMBLY.out.map { gca_value, last_id, extra -> tuple(gca_value, extra[0], last_id, extra[1]) }

    def update_keys_out = UPDATE_KEYS_METADATA(write2db_assembly_out)

    def write2db_metadata_in = update_keys_out.map { gca_value, metadata_json, species_tmp -> tuple(gca_value, metadata_json, [species_tmp]) }
    WRITE2DB_METADATA(write2db_metadata_in, write2db_script, false)
    def write2db_metadata_out = WRITE2DB_METADATA.out.map { gca_value, last_id, extra -> tuple(gca_value, extra[0], last_id) }

    def species_checker_out = SPECIES_CHECKER(write2db_metadata_out, species_checker_script)

    def write2db_species_in = species_checker_out.map { gca_value, species_json -> tuple(gca_value, species_json, []) }
    WRITE2DB_SPECIES(write2db_species_in, write2db_script, false)
    def write2db_species_out = WRITE2DB_SPECIES.out.map { gca_value, last_id, _extra -> tuple(gca_value, last_id) }

    def get_tolid_out = GET_TOLID(write2db_species_out)

    def write2db_tolid_in = get_tolid_out.map { gca_value, tolid_json -> tuple(gca_value, tolid_json, []) }
    WRITE2DB_TOLID(write2db_tolid_in, write2db_script, true)

    gca_list = WRITE2DB_TOLID.out.map { gca_value, _last_id, _extra -> gca_value }.map { it -> it.trim() }.collectFile(name: 'gca_list_to_report.txt', newLine: true)

    REPORT(gca_list, last_update)
}
