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

import java.nio.file.Files

// https://hub.docker.com/r/staphb/fastqc v12
process RUN_FASTQC {
    label 'fastqc'
    tag "$taxon_id:$run_accession"
    publishDir "${params.outDir}/$taxon_id/$run_accession/", mode: 'copy'
    afterScript "sleep $params.files_latency"
    input:
    tuple val(taxon_id), val(gca), val(run_accession), val(pair1), val(pair2), path(dataFileQuery)

    output:
    //val(fastqcOUT)
    tuple val(taxon_id), val(gca), val(run_accession), path(dataFileQuery),val("${params.outDir}/$taxon_id/$run_accession/fastqc")
    
    when:
    def file1 = new File(pair1)
    def file2 = new File(pair2)
    file1.exists() && file2.exists()

    script:
    """
    if [ ! -d "${params.outDir}/${taxon_id}/${run_accession}/fastqc/" ]; then
    mkdir -p fastqc 
    fastqc  ${pair1} ${pair2} --quiet --extract --threads ${task.cpus} --outdir fastqc; \
    cp -r fastqc ${params.outDir}/${taxon_id}/${run_accession}
    else
    echo "Directory fastqc already exists. Skipping process"
    fi
    """
}

