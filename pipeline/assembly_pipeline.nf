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

include { FETCH_GCA } from '../modules/fetch_gca.nf' 
include { PARSE_METADATA } from '../modules/parse_metadata.nf' 
include { WRITE2DB_METADATA } from '../modules/write2db_metadata.nf' 
include { UPDATE_KEYS_METRICS } from '../modules/update_keys_metrics.nf' 
include { WRITE2DB_METRICS } from '../modules/write2db_metrics.nf'
include { SPECIES_CHECKER } from '../modules/species_checker.nf' 
include { WRITE2DB_SPECIES } from '../modules/write2db_species.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW: REGISTER NEW ASSEMBLIES IN DB
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    FETCH_GCA(params.taxon, params.last_update)

    gca = FETCH_GCA.out.splitText().map{it -> it.trim()}

    PARSE_METADATA(gca)
    
    WRITE2DB_METADATA(PARSE_METADATA.out.gca, PARSE_METADATA.out.metadata, PARSE_METADATA.out.metrics_tmp, PARSE_METADATA.out.species_tmp)

    UPDATE_KEYS_METRICS(WRITE2DB_METADATA.out.gca, WRITE2DB_METADATA.out.metrics_tmp, WRITE2DB_METADATA.out.last_id, WRITE2DB_METADATA.out.species_tmp)

    WRITE2DB_METRICS(UPDATE_KEYS_METRICS.out.gca, UPDATE_KEYS_METRICS.out.metrics, UPDATE_KEYS_METRICS.out.species_tmp)

    SPECIES_CHECKER(WRITE2DB_METRICS.out.gca, WRITE2DB_METRICS.out.species_tmp)

    WRITE2DB_SPECIES(SPECIES_CHECKER.out.gca, SPECIES_CHECKER.out.species)

    /*
    WRITE2DB_SPECIES.out.gca.view()
    WRITE2DB_METRICS.out.last_id.view()
    */

}