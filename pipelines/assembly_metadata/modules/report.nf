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

process REPORT {
    publishDir "${params.output_dir}/nextflow_output/", mode: 'copy'

    input:
    path gca_list
    val last_update

    output:
    path "report.txt"
    path "gca_to_run_ncbi.csv"

    script:
    """
    chmod +x ../bin/python/create_report.py
    python ../bin/python/create_report.py \
    --file-list $gca_list --metadata ${params.metadata_params} --update-date $last_update
    """
}