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
    params.each{ k, v -> println "params.${k.padRight(25)} = ${v}" }



    FETCH_ASSEMBLIES(params.screen_date)
    def gca = FETCH_ASSEMBLIES.out.splitText().map{it -> it.trim()}

    FETCH_METADATA(gca)

    def attempt_update = FETCH_METADATA.out.attempt_update.map{it -> it.trim()}

    if (attempt_update) {
        def metadata_file = FETCH_METADATA.out.metadata_json..map{it -> it.trim()} //view { it -> println "json file: ${it}"}
    }

}

workflow.onComplete {
    log.info "Pipeline completed at: ${new Date().format('dd-MM-yyyy HH:mm:ss')}"
}

