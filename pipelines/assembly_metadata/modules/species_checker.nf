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
SPECIES CHECKER
This process checks the species information based on the provided JSON.
Inputs:
- gca: The GCA accession.
- species_tmp: The temporary species metadata file.
- last_id: The last processed ID.
Outputs:
- gca: The GCA accession.
- ${species_tmp.baseName}.json: The updated species metadata in JSON format.
*/

process SPECIES_CHECKER {

    label 'python'
    tag "$gca"
    publishDir "${params.output_dir}/nextflow_output/$gca", mode: 'copy'

    input:
    tuple val(gca), path(species_tmp), path(last_id)

    output:
    tuple val(gca), path("${species_tmp.baseName}.json")

    script:
    """
    species_checker.py --json-path $species_tmp --ncbi_url ${params.ncbi_url} --enscode ${params.enscode}
    """
}