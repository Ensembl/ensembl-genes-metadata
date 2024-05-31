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



// module description 
process RUN_STAR {
    tag "$run_accession"
    label 'star'
    publishDir "${params.outDir}/$taxon_id/$run_accession/star/", mode: 'copy'
    afterScript "sleep $params.files_latency"  // Needed because of file system latency

    input:
    tuple val(taxon_id), val(gca), val(run_accession), \
    val(fastqFile1), val(fastqFile2), val(genomeDir)

    output:
    tuple val(taxon_id), val(gca), val(run_accession),\
    val("${params.outDir}/${taxon_id}/${run_accession}/star/${run_accession}_Log.final.out")

    script:
    //def star_dir = new File(fastqFile1).parent
    def starTmpDir =  "${params.outDir}/${taxon_id}/${run_accession}/star/tmp"
    def outFileNamePrefix = "${params.outDir}/${taxon_id}/${run_accession}/star/${run_accession}_"
    """
    rm -rf ${starTmpDir}
    STAR --limitSjdbInsertNsj 2000000 \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --outSAMstrandField intronMotif --runThreadN ${task.cpus} \
    --twopassMode Basic --runMode alignReads --genomeDir ${file(genomeDir)} \
    --readFilesIn ${file(fastqFile1)} ${file(fastqFile2)} --outFileNamePrefix ${outFileNamePrefix} \
    --outSAMattrRGline "ID:${run_accession}" --outTmpDir ${starTmpDir} --outSAMtype BAM \
    SortedByCoordinate  --outBAMsortingBinsN 200 
    
    """
    //--readFilesCommand zcat  after subsampling they should be unzipped
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


/* INDEXING OPTION
def minSeqLength = 0
    def star_dir = "${params.outDir}/${taxon_id}/${run_accession}/star_alignment"
    
    def star_index_file = $star_dir + "/SAindex"
   
    # This calculates the base-2 logarithm of the genome_size. The logarithm of the genome size is
    # a measure of how many bits are needed to represent the genome size in binary.
    # The choice of 14 as the maximum value is likely based on empirical observations and optimization
    # considerations. Too large of a seed length can lead to increased memory usage and potentially
    # slower indexing, while a seed length that is too small might affect alignment accuracy.
     
    if !star_index_file.exists() :
        def seqRegionLength = getSeqRegionLength(genomeFile, minSeqLength)
        def genome_size = seqRegionToLength.values().sum()
        def index_bases = calculateIndexBases(${genome_size})  
        params.star_bin --runThreadN params.cpus --runMode "genomeGenerate" \
        --outFileNamePrefix $star_dir --genomeDir  $star_dir      \
        --genomeSAindexNbases $index_bases --genomeFastaFiles ${genome_file}

        
    def star_tmp_dir = star_dir / "tmp"

    params.star_bin --limitSjdbInsertNsj 2000000 \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --outSAMstrandField intronMotif --runThreadN params.cpus \
    --twopassMode Basic --runMode alignReads --genomeDir $genome_file \
    --readFilesIn $fastqFile1 $fastqFile2 --outFileNamePrefix $star_dir \
    --outTmpDir $tmp_dir --outSAMtype BAM SortedByCoordinate   

    """
}

import groovy.transform.TypeChecked

@TypeChecked
Map getSeqRegionLength(Path genomeFile, int minSeqLength) {
    def currentHeader = ""
    def currentSeq = ""
    def seqRegionToLength = [:]

    genomeFile.eachLine { line ->
        def match = line =~ />(.+)$/
        if (match && currentHeader) {
            if (currentSeq.size() > minSeqLength) {
                seqRegionToLength[currentHeader] = currentSeq.size()
            }
            currentSeq = ""
            currentHeader = match[0][1]
        } else if (match) {
            currentHeader = match[0][1]
        } else {
            currentSeq += line.trim()
        }
    }

    if (currentSeq.size() > minSeqLength) {
        seqRegionToLength[currentHeader] = currentSeq.size()
    }
    
    return seqRegionToLength
}

// Example usage:
// def genomeFile = new File("/path/to/your/genome_file.fa")
// def minSeqLength = 0
// def seqRegionLength = getSeqRegionLength(genomeFile, minSeqLength)
// println seqRegionLength
*/  
