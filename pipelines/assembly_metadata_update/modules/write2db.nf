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
WRITE2DB
Shared process that writes a taxonomy JSON file to the DB via write2db.py.
Used both for the integrity-check taxonomy update and the post taxonomy-check
update path.
Inputs:
- gca: The GCA accession.
- taxonomy_json: The taxonomy JSON file to write to the DB.
- passthrough: Opaque call-site data to forward unchanged alongside the output.
Outputs:
- gca: The GCA accession.
- ${taxonomy_json.baseName}.last_id: The last processed ID.
- passthrough: The call-site data, forwarded unchanged.
*/

process WRITE2DB {

    label 'python'
    tag "${gca}"
    publishDir "${params.output_dir}/nextflow_output/${gca}", mode: 'copy'

    input:
    tuple val(gca), path(taxonomy_json), val(passthrough)
    path write2db_script

    output:
    tuple val(gca), path("${taxonomy_json.baseName}.last_id"), val(passthrough)

    script:
    """
    python ${write2db_script} --file-path ${taxonomy_json} --metadata '${params.metadata_params_string}' --config ${params.db_table_conf}
    """
}
