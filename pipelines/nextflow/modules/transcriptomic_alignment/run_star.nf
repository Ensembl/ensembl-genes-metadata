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
    tuple val(taxon_id), val(gca), val(run_accession),\
    val("${params.outDir}/${taxon_id}/${run_accession}/star/${run_accession}_Log.final.out")

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
    --outSAMattrRGline "ID:${run_accession}" --outTmpDir ${starTmpDir} --outSAMtype BAM \
    SortedByCoordinate  
    """
    
}
/*
linuxbrew/bin/STAR  --limitBAMsortRAM 2006305390 --outBAMsortingThreadN 30
--limitSjdbInsertNsj 2000000 --outFilterIntronMotifs RemoveNoncanonicalUnannotated
--outSAMstrandField intronMotif --runThreadN 30 --twopassMode Basic 
--runMode alignReads --genomeDir /hps/nobackup/flicek/ensembl/genebuild/
ftricomi/fish/mummichog_annotation/fundulus_heteroclitus/GCA_011125445.2
/genome_dumps --readFilesIn /hps/nobackup/flicek/ensembl/genebuild/ftricomi
/fish/mummichog_annotation/fundulus_heteroclitus/GCA_011125445.2//rnaseq/input/
SRR12475462_1.fastq.gz /hps/nobackup/flicek/ensembl/genebuild/ftricomi/fish/
mummichog_annotation/fundulus_heteroclitus/GCA_011125445.2//rnaseq/input/S
RR12475462_2.fastq.gz --outFileNamePrefix /hps/nobackup/flicek/ensembl/genebuild/
ftricomi/fish/mummichog_annotation/fundulus_heteroclitus/GCA_011125445.2/rnaseq/
output/SRR12475462_   --outSAMattrRGline "ID:SRR12475462" 
--outTmpDir /hps/nobackup/flicek/ensembl/genebuild/ftricomi/fish/
mummichog_annotation/fundulus_heteroclitus/GCA_011125445.2/rnaseq/output
/SRR12475462_tmp --outSAMtype BAM SortedByCoordinate --outBAMsortingBinsN 200
*/
