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

process PROCESS_FASTQC_OUTPUT {
    label 'python'
    tag "$run_accession"
    publishDir "${params.outDir}/$taxon_id/$run_accession/fastqc", mode: 'copy'
    afterScript "sleep $params.files_latency"  // Needed because of file system latency
    
    input:
    tuple val(taxon_id), val(gca), val(run_accession), path(dataFileQuery),path(fastqc_dir), val(runId)

    output:
    tuple val(taxon_id), val(gca), val(run_accession)
    path("complete_insert_into_data_file.json")

    script:
    """
    chmod +x $projectDir/bin/parse_fastqc.py  # Set executable permissions
    parse_fastqc.py --fastqc_results_path ${fastqc_dir} --data_file_json ${dataFileQuery} --run_id ${runId}
    """
}

