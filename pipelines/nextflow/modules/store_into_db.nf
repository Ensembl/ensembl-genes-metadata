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
import java.nio.file.Files
import java.nio.file.Paths
import groovy.io.FileType
import java.util.concurrent.TimeUnit

include { setMetaDataRecord } from './utils.nf'
process STORE_INTO_DB {
    scratch false
    label 'default'
    tag "$run_accession"
    maxForks 10

    input:
    tuple val(taxon_id), val(gca), val(run_accession), val(query)
    
    output:
    tuple val(taxon_id), val(gca), val(run_accession)

    script:
    queriesArray = query.replaceAll('; ', ' ').split(";")
    log.info("Queries Array:  ${run_accession}  ${queriesArray} ")
    if (queriesArray.size() == 1) {
    setMetaDataRecord(queriesArray[0].trim().toString())
    } else if (queriesArray.size() > 1) {
    queriesArray.each { query ->
            // Trim the query to remove any leading/trailing whitespace
            query = query.trim()
            // Check if the query is not empty
            if (query) {
                log.info("STORE INTO DB  ${query};")
                setMetaDataRecord(query.toString())
            }
        }
    }

    """
    echo "${queriesArray.toString().trim()}"
    """       
}
