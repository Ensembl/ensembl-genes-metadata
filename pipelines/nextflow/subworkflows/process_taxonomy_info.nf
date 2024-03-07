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


includeConfig './nextflow.config'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/




def resultPresent = checkTableResult(jdbcUrl, username, password, query)
println "Result present: $resultPresent"

include { GET_RUN_ACCESSIONS } from '../modules/process_taxonomy_info/get_run_accessions.nf'


//  taxon id present or not? if yes get all new short read data after this date if not add it for the first time


workflow PROCESS_TAXONOMY_INFO {
    take:
    taxon_id                   

    // Define the output channel for run accessions
    output:
    val run_accessions_path into runAccessionsPath
    
    main:
    def taxonomyExists = checkTaxonomy(params.jdbcUrl, params.transcriptomic_user, params.transcriptomic_password, taxonId)
    if (taxonomyExists) {
      // Retrieve new run accessions for short-read transcriptomic data published AFTER the last check date
      def lastDate = getLastCheckedDate(params.jdbcUrl, params.transcriptomic_user, params.transcriptomic_password, taxonId)
      emit GET_RUN_ACCESSIONS (taxon_id, last_date)
      updateLastCheckedDate(params.jdbcUrl, params.transcriptomic_user, params.transcriptomic_password, taxonId)
    else{
      // Add the new taxon id and last_check=currentDate and retrieve all the run accessions for short-read transcriptomic data 
      def addTaxonId= insertMetaRecord(params.jdbcUrl, params.transcriptomic_user, params.transcriptomic_password, taxonId)
      emit GET_RUN_ACCESSIONS (taxon_id)
    }  
    
    }
}