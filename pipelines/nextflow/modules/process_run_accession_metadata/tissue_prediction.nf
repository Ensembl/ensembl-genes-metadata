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
process TISSUE_PREDICTION {
    label 'llm'
    tag "$taxon_id"
    //publishDir "${params.outDir}/$taxon_id/$run_accession", mode: 'copy'
    maxForks 1
    afterScript "sleep $params.files_latency"  // Needed because of file system latency
    
    input:
    tuple val(taxon_id)


    script:
    """
    # Check if Python dependencies are installed
    # Read each line in the requirements file
    while read -r package; do \\
    if ! pip show -q "\$package" &>/dev/null; then 
        echo "\$package is not installed" 
        pip install "\$package"
    else
        echo "\$package is already installed"
    fi
    done < ${projectDir}/bin/requirements_llm.txt

    chmod +x $projectDir/bin/llm_prediction.py  # Set executable permissions
    llm_prediction.py --taxon_id ${taxon_id}  --hugging_face_token  $params.hugging_face_token    --host $params.transcriptomic_dbhost --user $params.transcriptomic_dbuser --password $params.transcriptomic_dbpassword --database $params.transcriptomic_dbname --port $params.transcriptomic_dbport --batch_size $params.batch_size
    fi
    """
}


