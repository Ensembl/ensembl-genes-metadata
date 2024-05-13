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

include { checkTaxonomy } from '../utils.nf'
include { getLastCheckedDate } from '../utils.nf'
include { setMetaRecord } from '../utils.nf'

process PROCESS_TAXON_ID {

    label 'default'
    tag "$taxon_id"

    input:
    tuple val(taxon_id), val(gca)

    output:
    tuple(val(taxon_id), val(gca), stdout)


    script:
    //def lastCheckedDate =[]
    def taxonomyExists = checkTaxonomy(taxon_id)
    if (taxonomyExists){
        // Retrieve new run accessions for short-read transcriptomic data published AFTER the last check date
        lastCheckedDate = getLastCheckedDate(taxon_id)[0].last_check
        //updateLastCheckedDate(params.jdbcUrl, params.transcriptomic_dbuser, params.transcriptomic_dbpassword, taxonId)
    }else{
        // Add the new taxon id and last_check=currentDate and retrieve all the run accessions for short-read transcriptomic data 
        setMetaRecord(taxon_id,'insert')
        lastCheckedDate = '2019-01-01'
    }

    //lastCheckedDate.view()
    """
    echo $lastCheckedDate
    """
}    