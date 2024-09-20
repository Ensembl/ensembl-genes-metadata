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
process GET_RUN_ACCESSION_METADATA {
    label 'python'
    tag "$run_accession"
    publishDir "${params.outDir}/$taxon_id/$run_accession", mode: 'copy'
    maxForks 25
    afterScript "sleep $params.files_latency"  // Needed because of file system latency
    
    input:
    tuple val(taxon_id), val(gca), val(run_accession)

    output:
    tuple val(taxon_id), val(gca), val(run_accession), path("insert_into_run.json")
    tuple val(taxon_id), val(gca), val(run_accession), path("insert_into_study.json")
    tuple val(taxon_id), val(gca), val(run_accession), path("insert_into_data_file.json")

    script:
    """
    if [ ! -s "${params.outDir}/${taxon_id}/${run_accession}/insert_into_run.json" ] || \
        [ ! -s "${params.outDir}/${taxon_id}/${run_accession}/insert_into_study.json" ] || \
        [ ! -s "${params.outDir}/${taxon_id}/${run_accession}/insert_into_data_file.json" ]; then

    # Check if Python dependencies are installed
    # Read each line in the requirements file
    while read -r package; do \\
    if ! pip show -q "\$package" &>/dev/null; then 
        echo "\$package is not installed" 
        pip install "\$package"
    else
        echo "\$package is already installed"
    fi
    done < ${projectDir}/bin/requirements.txt

    chmod +x $projectDir/bin/get_metadata.py  # Set executable permissions
    get_metadata.py --run ${run_accession}
    else
        echo "Metadata already exist.";
        cp -r ${params.outDir}/${taxon_id}/${run_accession}/insert_into* .
    fi
    """
}
