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
TAXONOMY CHECK
This process checks the taxonomy information for a given GCA accession using the taxonomy.py script.
Inputs:
- gca: The GCA accession for which to check taxonomy information.
- attempt_update: A flag indicating whether to attempt an update (true/false).
- metadata_json: The JSON file containing the metadata information.
Outputs:
- OLD_TAXON_ID: The previous taxon ID before the update.
- NEW_TAXON_ID: The new taxon ID after the update.
- STATUS: The status of the taxonomy check.
*/

process TAXONOMY_CHECK {

    label 'python'
    tag "${gca}"

    input:
    tuple val(gca), val(attempt_update), path(metadata_json)

    output:
    tuple val(gca), val(attempt_update), path(metadata_json), env('OLD_TAXON_ID'), env('NEW_TAXON_ID'), env('STATUS'), emit: taxonomy_check

    when:
    attempt_update.trim() == 'true'

    script:
    """
    # Run the script and capture output
    OUTPUT=\$(taxonomy.py --accession_json ${metadata_json} --accession ${gca} \
             --metadata_params '${params.metadata_params_string}' --taxonomy_check)

    # Parse the comma-separated output (old_taxon_id,new_taxon_id,status)
    IFS=',' read -r OLD_TAXON_ID NEW_TAXON_ID STATUS <<< "\$OUTPUT"
    export OLD_TAXON_ID NEW_TAXON_ID STATUS
    """
}
