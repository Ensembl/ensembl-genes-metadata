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
ASSEMBLY NAME
This process runs the assembly_name.py script to update assembly names in the database. 
Inputs:
- gca: The GCA accession number for the assembly.
- attempt_update: A flag indicating whether to attempt the update (true/false).
- metadata_json: Path to the JSON file containing metadata for the assembly.
Outputs:
- asm_name_update: The standard output from the assembly_name.py script, a string indicating the result of the update operation.
*/

process ASSEMBLY_NAME {
    
    label 'python'
    tag "$gca"

    input:
    tuple val(gca), val(attempt_update), path(metadata_json)

    output:
    stdout emit: asm_name_update

    when:
    attempt_update.trim() == 'true'

    script:
    """
    assembly_name.py --accession_json $metadata_json --accession $gca  --metadata_params '${params.metadata_params_string}'
    """
}
