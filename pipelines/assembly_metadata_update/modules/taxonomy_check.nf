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


process TAXONOMY_CHECK {
    
    label 'python'
    tag "$gca"
    
    when:
    attempt_update.trim() == 'true'

    input:
    tuple val(gca), val(attempt_update), path(metadata_json)

    output:
    tuple val(gca), val(attempt_update), path(metadata_json), 
          env(OLD_TAXON_ID), env(NEW_TAXON_ID), env(STATUS), 
          emit: taxonomy_check

    script:
    """
    # Run the script and capture output
    OUTPUT=\$(taxonomy.py --accession_json $metadata_json --accession $gca \
             --metadata_params ${params.metadata_params} --taxonomy_check)
    
    # Parse the comma-separated output
    export OLD_TAXON_ID=\$(echo \$OUTPUT | cut -d',' -f1)
    export NEW_TAXON_ID=\$(echo \$OUTPUT | cut -d',' -f2)
    export STATUS=\$(echo \$OUTPUT | cut -d',' -f3)
    """
}

// process TAXONOMY_CHECK {
    
//     label 'python'
//     tag "$gca"
    
//     when:
//     attempt_update.trim() == 'true'

//     input:
//     tuple val(gca), val(attempt_update), path(metadata_json)

//     output:
//     tuple val(gca), val(attempt_update), path(metadata_json), stdout, emit: taxonomy_check

//     script:
//     """
//     taxonomy.py --accession_json $metadata_json --accession $gca  --metadata_params ${params.metadata_params} --taxonomy_check
//     """
// }
