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

include { GET_RUN_ACCESSION } from '../modules/process_taxonomy_info/check_run_accession.nf'
include { CHECK_RUN_ACCESSION } from '../modules/process_taxonomy_info/get_run_accession.nf'



workflow PROCESS_TAXONOMY_INFO {
    take:
    taxon_id    
    run_accession            
    transcriptomic_dbname                 
    transcriptomic_host                  
    transcriptomic_port   
    transcriptomic_user    
    transcriptomic_password           

    main:

    filtered_run_accessions = Channel.empty()
    inital_run_accession = GET_RUN_ACCESSION (taxon_id)  ///list of original run accession for a taxon id
    //now we need to filter them
    filtered_run_accessions = CHECK_RUN_ACCESSION (
      taxon_id,
      transcriptomic_dbname, 
      transcriptomic_host,
      transcriptomic_port,   
      transcriptomic_user,
      transcriptomic_password,
      inital_run_accession.out.run_accession_list)


    emit:
    filtered_run_accessions            = filtered_run_accessions                  // channel: [run_accessions]
}