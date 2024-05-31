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
    RUN STAR ALIGNER
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FETCH_GENOME } from '../modules/transcriptomic_alignment/fetch_genome.nf'
include { INDEX_GENOME } from '../modules/transcriptomic_alignment/index_genome.nf'
include { RUN_STAR } from '../modules/transcriptomic_alignment/run_star.nf'
include { EXTRACT_UNIQUELY_MAPPED_READS_PERCENTAGE } from '../modules/transcriptomic_alignment/extract_uniquely_mapped_reads_percentage.nf'
include { BUILD_QUERY as BUILD_QUERY_ALIGN_METADATA } from '../modules/build_query.nf'
include { CLEANING } from '../modules/cleaning.nf'

include { setMetaDataRecord } from '../modules/utils.nf'
include { getDataFromTable } from '../modules/utils.nf'
include { updateTable } from '../modules/utils.nf'


workflow RUN_ALIGNMENT{
    take:
    shortReadMetadata

    main:
    def genomeAndDataToAlign = FETCH_GENOME(shortReadMetadata)
    def genomeIndexAndDataToAlign = INDEX_GENOME(genomeAndDataToAlign)
    def starOutput = RUN_STAR(genomeIndexAndDataToAlign)
    def (starMetadata, insertIntoAlign) = EXTRACT_UNIQUELY_MAPPED_READS_PERCENTAGE(starOutput)
    def updateValue = "False"
    def (runAccessionData_StarOutput,insertIntoAlignQuery) = BUILD_QUERY_ALIGN_METADATA(starMetadata, insertIntoAlign, updateValue)
    def data3=insertIntoAlignQuery
    data3.flatten().view { d -> "query ${d}"}
    insertIntoAlignQuery.subscribe { line ->
        setMetaDataRecord(line.toString())
    }
    def runAccessionData_NewQCstatus = runAccessionData_StarOutput.map { result ->
        def (taxon_id, gca, run_accession) = result
        def run_Id = getDataFromTable("run_id", "run", "run_accession", run_accession)[0].run_id
        //updateRunStatus(runId)
        updateTable("run_id", run_Id, "run", "qc_status", "ALIGNED")
        setLastCheckDate(taxon_id,'update')
        return tuple(taxon_id, gca, run_accession)
    }
    CLEANING(runAccessionData_NewQCstatus)

}

