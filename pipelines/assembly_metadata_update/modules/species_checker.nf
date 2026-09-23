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
This process checks the species information for a given GCA accession using the species_checker.py script.
Inputs:
- gca: The GCA accession for which to check species information.
- attempt_update: A flag indicating whether to attempt an update (true/false).
- metadata_json: The JSON file containing the metadata information.
- old_taxon_id: The previous taxon ID before the update.
- new_taxon_id: The new taxon ID after the update.
Outputs:
- taxonomy_${new_taxon_id}.json: The JSON file containing the updated taxonomy information.
*/

process SPECIES_CHECKER {

    label 'python'
    tag "${gca}"
    publishDir "${params.output_dir}/nextflow_output/${gca}", mode: 'copy'

    input:
    tuple val(gca), val(attempt_update), path(metadata_json), val(old_taxon_id), val(new_taxon_id)
    path species_checker_script

    output:
    tuple val(gca), val(attempt_update), path(metadata_json), path("taxonomy_${new_taxon_id}.json")

    when:
    attempt_update.trim() == 'true'

    script:
    """
    python ${species_checker_script} --taxon_id ${new_taxon_id} --ncbi_url ${params.ncbi_url} --taxonomy_update
    """
}
