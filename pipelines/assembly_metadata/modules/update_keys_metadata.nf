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
UPDATE_KEYS_METADATA
This process updates the keys in the metadata JSON based on the provided temporary metadata and last processed ID
Inputs:
- gca: The GCA accession.
- metadata_tmp: The temporary metadata file.
- last_id: The last processed ID.
- species_tmp: The temporary species metadata file.
Outputs:
- gca: The GCA accession.
- ${metadata_tmp.baseName}.json: The updated metadata in JSON format.
- species_tmp: The temporary species metadata file (unchanged).
*/

process UPDATE_KEYS_METADATA {

    label 'python'
    tag "$gca"
    publishDir "${params.output_dir}/nextflow_output/$gca", mode: 'copy'
    
    input:
    tuple val(gca), path(metadata_tmp), path(last_id), path(species_tmp)
    
    output:
    tuple val(gca), path("${metadata_tmp.baseName}.json"), path(species_tmp)
    
    script:
    """
    update_keys.py --json-path $metadata_tmp --file-id-path $last_id --config ${params.db_table_conf}
    """
}
