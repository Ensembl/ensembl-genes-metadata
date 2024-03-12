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

// https://hub.docker.com/r/staphb/fastqc v12
process RUN_FASTQC {
    scratch true
    label 'fastqc'
    tag "$taxon_id"

    input:
    val taxon_id
    val run accession
    set pair1, pair2 from pairedFastqFiles

    output:
    storeDir joinPath(params.outDir, "${taxon_id}", "${run_accession}", "fastqc") into fastqcOutput


    script:
    """
    fastqc fastqc_output ${pair1} ${pair2} --quiet --extract --threads ${task.cpus}
    """

    emit:
    emit fastqcOutput
}

