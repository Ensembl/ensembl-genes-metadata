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

process SAM2BAM {
    tag "$run_accession"
    label 'samtools'
    storeDir "${params.outDir}/$taxon_id/$run_accession/alignment/"
    afterScript "sleep $params.files_latency"  // Needed because of file system latency

    input:
    tuple val(taxon_id), val(genomeDir), val(tissue), val(run_accession), path(sam_file)


    output:
    //tuple val(taxon_id), val(genomeDir), val(gca), val(platform), val(paired), val(tissue), val(run_accession), path("*.bam")
    tuple val(taxon_id), val(genomeDir), val(tissue), val(run_accession), path("*.bam")




    script:
    """
    samtools index ${sam_file} ${run_accession}.bam
    samtools index ${run_accession}.bam ${run_accession}.bam.bai
    """
}


