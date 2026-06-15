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
/*
PARSE_METADATA
This process retrieves the assembly metadata for a given GCA accession from the NCBI
Inputs:
- gca: The GCA accession for which to retrieve assembly metadata.
Outputs:
- gca: The GCA accession.
- ${gca}_assembly.json: The assembly metadata in JSON format.
- ${gca}_metadata.tmp: The assembly metadata in temporary format.
- ${gca}_species.tmp: The species metadata in temporary format. 
*/

process PARSE_METADATA {
    
    label 'python'
    tag "$gca"
    publishDir "${params.output_dir}/nextflow_output/$gca", mode: 'copy'

    input:
    val gca

    output:
    tuple val(gca), path("${gca}_assembly.json"), path("${gca}_metadata.tmp"), path("${gca}_species.tmp")

    script:
    """
    retrieving_metadata.py --accession $gca --ncbi_url ${params.ncbi_url}
    """
}
