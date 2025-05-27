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




process DOWNLOAD_FASTQS {
    label "python"
    tag "${taxon_id}:${run_accession}"
    maxForks 25
    storeDir "${params.outDir}/$taxon_id/$run_accession"
    afterScript "sleep $params.files_latency"  // Needed because of file system latency
    conda "${projectDir}/bin/environment.yml"

    input:
    tuple val(taxon_id), val(gca), val(platform), val(paired), val(tissue), val(run_accession), val(genomeDir),  val(url1), val(md5_1), val(url2),  val(md5_2)

    output:
    tuple val(taxon_id), val(genomeDir), val(gca), val(platform), val(paired), val(tissue), val(run_accession), path("${params.outDir}/$taxon_id/$run_accession/*_1.fastq.gz") , val(paired ? path("${params.outDir}/$taxon_id/$run_accession/*_2.fastq.gz") : null)
//input:
//    tuple val(taxon_id), val(gca), val(run_accession), path(pair1), path(pair2, optional: true), val(genomeDir)

    script: 
    """"
    chmod +x $projectDir/bin/download_fastqs.py
    python3 download_fastqs.py \\
        --taxon_id ${taxon_id} \\
        --run_accession ${run_accession} \\
        --url1 ${url1} \\
        --md5_1 ${md5_1} \\
        ${paired ? "--url2 ${url2} --md5_2 ${md5_2} --paired" : ""} \\
        --outDir ${params.outDir} \\

    """
    
}
