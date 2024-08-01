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


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BUILD_QUERY as BUILD_QUERY_RUN_METADATA } from '../modules/build_query.nf'
include { BUILD_QUERY as BUILD_QUERY_STUDY_METADATA} from '../modules/build_query.nf'
include { BUILD_QUERY as BUILD_QUERY_DATA_FILE_METADATA} from '../modules/build_query.nf'
include { STORE_INTO_DB as STORE_INTO_DB_RUN} from '../modules/store_into_db.nf'
include { STORE_INTO_DB as STORE_INTO_DB_STUDY} from '../modules/store_into_db.nf'
include { STORE_INTO_DB as STORE_INTO_DB_DATA_FILE } from '../modules/store_into_db.nf'
include { DOWNLOAD_PAIRED_FASTQS } from '../modules/process_run_accession_metadata/download_paired_fastqs.nf'
include { GET_RUN_ACCESSION_METADATA } from '../modules/process_run_accession_metadata/get_run_accession_metadata.nf'

include { PROCESS_FASTQC_OUTPUT } from '../modules/fastqc_processing/process_fastqc_output.nf'
include { RUN_FASTQC } from '../modules/fastqc_processing/run_fastqc.nf'
include { BUILD_QUERY } from '../modules/build_query.nf'
include { STORE_INTO_DB } from '../modules/store_into_db.nf'
include { SUBSAMPLE_FASTQ_FILES } from '../modules/fastqc_processing/subsample_fastq_files.nf'
include { ADAPTER_TRIMMING } from '../modules/fastqc_processing/adapter_trimming.nf'
include { CLEANING } from '../modules/cleaning.nf'

include { setMetaDataRecord } from '../modules/utils.nf'
include { checkRunFastQCStatus } from '../modules/utils.nf'
include { checkOverrepresentedSequences } from '../modules/utils.nf'
include { getDataFromTable } from '../modules/utils.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: PROCESS RUN ACCESSION METADATA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PROCESS_RUN_ACCESSION_METADATA {
    take:
    run_accession_list
    //tuple val(taxon_id), val(gca), val(run_accession)
    main:
    def db_meta1=run_accession_list
    db_meta1.flatten().view { d -> "GCA: ${d.gca}, Taxon ID: ${d.taxon_id}, run: ${d.run_accession}"}
    def(insertIntoRun, insertIntoStudy, queryDataFile) = GET_RUN_ACCESSION_METADATA(run_accession_list.flatten())
    def updateValue = "False"
    //insert in study table
    def (runAccessionMetadataStudyQuery) = BUILD_QUERY_STUDY_METADATA(insertIntoStudy, updateValue)
    def (runAccessionMetadataStudyOutput) = STORE_INTO_DB_STUDY(runAccessionMetadataStudyQuery)

    def (runAccessionMetadataRunQuery) = BUILD_QUERY_RUN_METADATA(insertIntoRun, updateValue)
    //STORE_INTO_DB_RUN(runAccessionMedatadataRunQuery)
    def (runAccessionMetadataRunOutput) = STORE_INTO_DB_RUN(runAccessionMetadataRunQuery)
    //insert in run table
    //insertIntoRunQuery.subscribe { line ->
    //    setMetaDataRecord(line.toString())
    //}
    //insert in study table
    //def (runAccessionMedatadataStudyQuery) = BUILD_QUERY_STUDY_METADATA(insertIntoStudy, updateValue)
    //def (runAccessionMedatadataStudyOutput) = STORE_INTO_DB_STUDY(runAccessionMetadataStudyQuery)
    //insertIntoStudyQuery.subscribe { line ->
    //    setMetaDataRecord(line.toString())
    //}
    def pairedFastqFiles=DOWNLOAD_PAIRED_FASTQS(runAccessionMetadataRunOutput)
    //def pairedFastqFilesMetadata = pairedFastqFiles.map { result ->
    //    def (taxon_id, gca, run_accession,pair1,pair2,dataFileQuery) = result
    //    log.info("resultsDOWNLOAD: ${result} ${result.gca[0]}")
    //    return tuple(result.taxon_id[0], result.gca[0], result.run_accession[0], result.pair1[0], result.pair2[0], result.dataFileQuery[0])
    //    }
    def fastqcOutput = RUN_FASTQC(pairedFastqFiles)
    def processedFastQCOutput = fastqcOutput.map { result ->
        def (taxon_id, gca, run_accession, dataFileQuery, fastqcPath) = result
        log.info("results: ${result}")
        log.info("Run Accession: ${run_accession}")
        //log.info("Run Accession: ${result.run_accession[0].toString()}")
        def run_Id = getDataFromTable("run_id", "run", "run_accession",  run_accession)[0].run_id
        //def run_Id = getDataFromTable("run_id", "run", "run_accession",  result.run_accession[0].toString())[0].run_id
        log.info "Run ID: ${run_Id}"
        return tuple(taxon_id, gca, run_accession, dataFileQuery, fastqcPath, run_Id.toString())

   //     return tuple(result.taxon_id[0], result.gca[0], result.run_accession[0], result.dataFileQuery[0], result.fastqcDir[0], run_Id.toString())
     }
    //def cc=processedFastQCOutput
    //cc.flatten().view { d -> "GCAFASTQC: ${d.gca}, Taxon ID: ${d.taxon_id}, run: ${d.run_accession}"}
     def (insertIntoDataFile) = PROCESS_FASTQC_OUTPUT(processedFastQCOutput)
     def (runAccessionDataFileQuery) = BUILD_QUERY_DATA_FILE_METADATA(insertIntoDataFile, updateValue)
     def (runAccessionDataFileOutput) = STORE_INTO_DB_DATA_FILE(runAccessionDataFileQuery)
    //insertIntoDataFileQuery.subscribe { line ->
     //    log.info("insertIntoDataFileQuery.subscribe ${line}")
     //    def queriesArray = line.toString().split(";")
     //    queriesArray.eachWithIndex { query, index ->
         // Trim the query to remove any leading/trailing whitespace
     //    query = query.trim()
         // Check if the query is not empty
      //   if (query) {
     //        log.info("queriesArray[${index}] ${query};")
      //       setMetaDataRecord(query.toString() + ";")                                                                   //  }
                                                                                                            //            }
                                                                                                          //           }
     def runAccessionData_QCstatus = runAccessionDataFileOutput.map { result ->
         def (taxon_id, gca, run_accession) = result
         log.info("runAccessionData_QCstatus ${run_accession}")
         def run_Id = getDataFromTable("run_id", "run", "run_accession", run_accession.toString())[0].run_id.toString()
             checkRunFastQCStatus(run_Id)
             return tuple(taxon_id, gca, run_accession)
         }
         //if (qc_status == 'QC_PASS') {
         def subsamplingOutput = SUBSAMPLE_FASTQ_FILES(runAccessionData_QCstatus)
         def subsampling_Output = subsamplingOutput.map { subsampling ->
         def (taxon_id, gca, run_accession, sub1, sub2) = subsampling
         log.info("subsampling: ${subsampling}")
         return [taxon_id:taxon_id, gca:gca, run_accession:run_accession, pair1:sub1, pair2:sub2]
     }

    
    
    emit:
    //pairedFastqFiles_metadata =   pairedFastqFiles      // channel: [taxon_id, gca, run_accession, path1, path2, data_file_query]
    subsamplingOutputMetadata = subsampling_Output //channel: [taxon_id, gca, run_accession, path1, path2]
}
