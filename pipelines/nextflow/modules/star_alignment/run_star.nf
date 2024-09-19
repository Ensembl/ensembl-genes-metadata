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

process RUN_STAR {
    tag "$run_accession"
    label 'star'
    publishDir "${params.outDir}/$taxon_id/$run_accession/star/", mode: 'copy'
    afterScript "sleep $params.files_latency"  // Needed because of file system latency

    input:
    tuple val(taxon_id), val(gca), val(run_accession), \
    val(pair1), val(pair2), val(genomeDir)

    output:
    tuple val(taxon_id), val(gca), val(run_accession), val(genomeDir),\
    val("${params.outDir}/${taxon_id}/${run_accession}/star/${run_accession}_Aligned.sortedByCoord.out.bam")
    // val("${params.outDir}/${taxon_id}/${run_accession}/star/${run_accession}_SJ.out.tab")

    script:
    def starTmpDir =  "${params.outDir}/${taxon_id}/${run_accession}/star/tmp"
    def outFileNamePrefix = "${params.outDir}/${taxon_id}/${run_accession}/star/${run_accession}_"
    """
    rm -rf ${starTmpDir}
    STAR --limitSjdbInsertNsj 2000000 \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --outSAMstrandField intronMotif --runThreadN ${task.cpus} \
    --twopassMode Basic --runMode alignReads --genomeDir ${file(genomeDir)} \
    --readFilesIn ${pair1} ${pair2} --outFileNamePrefix ${outFileNamePrefix} \
    --readFilesCommand zcat --outSAMattrRGline "ID:${run_accession}" --outTmpDir ${starTmpDir} --outSAMtype BAM \
    SortedByCoordinate  
    """
    
}
