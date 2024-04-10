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

includeConfig './pipelines/workflows/nextflow.config'
include { calculateIndexBases} from './pipelines/modules/utils.nf'


// module description 
process RUN_STAR {
    scratch false
    label 'run_star'
    tag 'star'
    
    input:
    tuple val(gca),  file(genome_file), file(fastqFile1), file(fastqFile2)

    output:
    file "${params.outDir}/${taxon_id}/${run_accession}/star_alignment/*.fastq.gz.Log.final.out", emit:log_final_out

    script:
    """
    def minSeqLength = 0
    def star_dir = "${params.outDir}/${taxon_id}/${run_accession}/star_alignment"
    
    def star_index_file = $star_dir + "/SAindex"
    /*
    # This calculates the base-2 logarithm of the genome_size. The logarithm of the genome size is
    # a measure of how many bits are needed to represent the genome size in binary.
    # The choice of 14 as the maximum value is likely based on empirical observations and optimization
    # considerations. Too large of a seed length can lead to increased memory usage and potentially
    # slower indexing, while a seed length that is too small might affect alignment accuracy.
    */   
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
