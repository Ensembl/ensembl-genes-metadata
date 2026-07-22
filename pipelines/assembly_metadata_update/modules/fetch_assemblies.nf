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
FETCH ASSEMBLIES
This process receive a list of GCAS accessions from the user or fetches a list from the assembly metadata database based on a screem date.
Inputs:
- screen_date: The date to filter assemblies in the database.
- gca_input: A flag indicating whether to use a user-provided list of GCA accessions (true/false).
- gca_list: Path to the file containing a list of GCA accessions (if gca_input is true).
Outputs:
- stdout: The standard output from the fetch_assemblies.py script or the grep command, which is a list of GCA accessions to be processed.
*/

process FETCH_ASSEMBLIES {

    label 'python'
    tag "date:$screen_date"

    input:
    val screen_date

    output:
    stdout

    script:
    if (params.gca_input){
    """
        grep '^GCA_' ${params.gca_list} 
    """
    }
    else {
    """
    fetch_assemblies.py --metadata ${params.metadata_params} --screen_date $screen_date
    """
    }

}