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
Shared process that writes a single JSON/tmp file to the DB via write2db.py.
Used at each stage of the registration pipeline (assembly, metadata, species,
tolid); files that need to ride alongside for the next process in the chain
are threaded through as an opaque passthrough path list rather than being
known to this process.
Inputs:
- gca: The GCA accession.
- file_to_write: The file to write to the DB.
- passthrough: Other files (possibly empty) to forward unchanged alongside the output.
- update_flag: Whether to pass --update to write2db.py.
Outputs:
- gca: The GCA accession.
- ${file_to_write.baseName}.last_id: The last processed ID.
- passthrough: The other files, forwarded unchanged.
*/

process WRITE2DB {

    label 'python'
    tag "${gca}"
    publishDir "${params.output_dir}/nextflow_output/${gca}", mode: 'copy'

    input:
    tuple val(gca), path(file_to_write), val(passthrough)
    path write2db_script
    val update_flag

    output:
    tuple val(gca), path("${file_to_write.baseName}.last_id"), val(passthrough)

    script:
    def flag = update_flag ? '--update' : ''
    """
    python ${write2db_script} --file-path ${file_to_write} ${flag} --metadata ${params.metadata_params} --config ${params.db_table_conf}
    """
}
