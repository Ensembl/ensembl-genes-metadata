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
REPORT UPDATE
This process updates the report for a given GCA accession using the report_update.py script.
Inputs:
- gca: The GCA accession for which to update the report.
- check: The type of check being performed.
- reporting: A flag indicating whether reporting is enabled (true/false).
- old_value: The previous value before the update.
- new_value: The new value after the update.
Outputs:
- slack_reporting.log: The log file containing the Slack reporting information.
*/

process REPORT_UPDATE {

    label 'python'
    tag "${gca}"
    publishDir "${params.output_dir}/nextflow_output/${gca}", mode: 'copy'

    input:
    tuple val(gca), val(check), val(reporting), val(old_value), val(new_value)

    output:
    path "slack_reporting.log", emit: asm_file

    when:
    reporting.trim() == 'true'

    script:
    """
    report_update.py \
    --accession ${gca} \
    --check_type ${check} \
    --previous_value "${old_value}" \
    --current_value "${new_value}" \
    --slack_users ${params.slack_user} \
    --metadata_params '${params.metadata_params_string}' \
    --slack_params '${params.slack_params}'
    """
}
