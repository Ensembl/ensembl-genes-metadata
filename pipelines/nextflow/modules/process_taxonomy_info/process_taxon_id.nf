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

include { setLastCheckDate } from '../utils.nf'
include { getDataFromTable } from '../utils.nf'

process PROCESS_TAXON_ID {

    label 'default'
    tag "$taxon_id"

    input:
    tuple val(taxon_id), val(gca)

    output:
    tuple(val(taxon_id), val(gca), stdout)

    script:
    def taxonomyExists = getDataFromTable('taxon_id', 'meta', 'taxon_id', taxon_id)
    if (!taxonomyExists.isEmpty()) {  
        lastCheckedDate = getDataFromTable('last_check', 'meta',  'taxon_id',  taxon_id)[0].last_check
    }else{
        // Add the new taxon id and last_check='2019-01-01'
        setLastCheckDate(taxon_id,'insert')
        lastCheckedDate = '2019-01-01'
    }

    """
    echo $lastCheckedDate
    """
}    
