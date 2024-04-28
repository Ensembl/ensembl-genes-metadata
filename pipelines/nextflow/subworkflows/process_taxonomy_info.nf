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
include { checkTaxonomy, getLastCheckedDate, updateLastCheckedDate, insertMetaRecord } from './pipelines/modules/utils.nf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GET_RUN_ACCESSIONS } from '../modules/process_taxonomy_info/get_run_accessions.nf'
include { PROCESS_TAXON_ID } from '../modules/process_taxonomy_info/process_taxon_id.nf'

workflow PROCESS_TAXONOMY_INFO {
    take:
    val data                   

    // Define the output channel for run accessions
    output:
    val(runAccessionList)
    
    main:
    def data1=data
    data1.flatten().view { d -> "Taxon ID: ${d.taxon_id}, GCA: ${d.gca},"}
    def taxonomyInfo = PROCESS_TAXON_ID(data.flatten())
    def (runAccessionList, runAccessionFile) = GET_RUN_ACCESSIONS(taxonomyInfo)
    
}