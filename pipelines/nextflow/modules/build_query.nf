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


process BUILD_QUERY {
    scratch false
    label 'python'
    tag "$run_accession"
    
    input:
    tuple val(taxon_id), val(gca), val(run_accession)
    path metadata2process
    val mysqlUpdate
    
    output:
    tuple val(taxon_id), val(gca), val(run_accession)
    stdout

    script:
    log.info("Executing Python script to build query for run: $run_accession")
    """
    chmod +x $projectDir/bin/write2db.py  # Set executable permissions
    write2db.py --file-path ${metadata2process} --update ${mysqlUpdate}
    """
}
