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
TAXONOMY UPDATE
This process updates the taxonomy information for a given GCA accession using the taxonomy.py script.
Inputs:
- gca: The GCA accession for which to update taxonomy information.
- attempt_update: A flag indicating whether to attempt an update (true/false).
- metadata_json: The JSON file containing the metadata information.
Outputs:
- stdout: The standard output from the taxonomy.py script.
*/

process TAXONOMY {

    label 'python'
    tag "${gca}"

    input:
    tuple val(gca), val(attempt_update), path(metadata_json)

    output:
    stdout

    when:
    attempt_update.trim() == 'true'

    script:
    """
    taxonomy.py --accession_json ${metadata_json} --accession ${gca}  --metadata_params '${params.metadata_params_string}' --taxonomy_update
    """
}
