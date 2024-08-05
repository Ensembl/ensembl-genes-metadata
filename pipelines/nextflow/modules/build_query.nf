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
include { setMetaDataRecord } from './utils.nf'
process BUILD_QUERY {
    scratch false
    label 'python'
    tag "$run_accession"
    //maxForks 50

    input:
    tuple val(taxon_id), val(gca), val(run_accession), path(metadata2process)
    val mysqlUpdate
    
    output:
    tuple val(taxon_id), val(gca), val(run_accession), stdout

    script:
    """
    chmod +x $projectDir/bin/write2db.py;
    $projectDir/bin/write2db.py --file_path ${metadata2process} --update ${mysqlUpdate} 
    """            
}
