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

process REPORT_UPDATE {
    
    label 'python'
    tag "$gca"
    publishDir "${params.output_dir}/nextflow_output/$gca", mode: 'copy'

    input:
    tuple val(gca), val(check), val(reporting), val(old_value), val(new_value)
    
    output:
    path "slack_reporting.log", emit: asm_file

    when:
    reporting.trim() == 'true'
    
    script:
    """
    report_update.py \
    --accession $gca \
    --check_type $check \
    --previous_value "$old_value" \
    --current_value "$new_value" \
    --slack_users $params.slack_user \
    --metadata_params '$params.metadata_params_string' \
    --slack_params '$params.slack_params'
    """
}
