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
INTEGRITY CHECKER
This process checks the integrity of the metadata for a given GCA accession using the integrity_checker.py script.
Inputs:
- gca: The GCA accession for which to check integrity.
Outputs:
- stdout: The standard output from the integrity_checker.py script. 
*/

process INTEGRITY_CHECKER {

    label 'python'
    tag "${gca}"

    input:
    val gca

    output:
    tuple val(gca), stdout

    script:
    """
    integrity_checker.py  --accession ${gca} --metadata '${params.metadata_params_string}' --delete
    """
}
