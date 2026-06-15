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
GET_TOLID
This process retrieves the TOLID for a given GCA accession from id.tol.sanger.ac.uk
Inputs:
- gca: The GCA accession for which to fetch the TOLID.
Outputs:
- gca_tolid.json: A JSON file containing the TOLID information for the given GCA accession.
*/

process GET_TOLID {

    label 'python'
    tag "$gca"
    publishDir "${params.output_dir}/nextflow_output/$gca", mode: 'copy'

    input:
    tuple val(gca), path(last_id)

    output:
    tuple val(gca), path("${gca}_tolid.json")

    script:
    """
    get_tolid.py --accession $gca --metadata ${params.metadata_params}
    """


}