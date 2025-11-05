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

process FETCH_GCA {

    label 'python'
    tag "taxon:$taxon"

    input:
    val taxon
    val last_update

    output:
    stdout

    script:
    if (params.add_gca) {
    """
        grep '^GCA_' ${params.gca_list}
    """
    }
    else {
    """
    fetch_new_assemblies.py \
    --taxon $taxon --date_update $last_update --db asm_metadata \
    --registry ${params.registry_params} --metadata ${params.metadata_params} \
    --ncbi ${params.ncbi_params} --ncbi_url ${params.ncbi_url}
    """
    }

}