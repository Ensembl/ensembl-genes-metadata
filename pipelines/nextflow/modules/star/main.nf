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

/* STAR documentation https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf
*/

process STAR {
    tag "$run_accession"
    label 'star'
    publishDir "${params.outDir}/$taxon_id/$run_accession/alignment/", mode: 'copy'
    afterScript "sleep $params.files_latency"  // Needed because of file system latency

    input:
    tuple val(taxon_id), val(genomeDir), val(platform), val(tissue), val(run_accession), val(pair1), val(pair2)
    output:
    //tuple val(taxon_id), val(genomeDir), val(gca), val(platform), val(paired), val(tissue), val(run_accession), path("*_Aligned.sortedByCoord.out.bam")
    tuple val(taxon_id), val(genomeDir), val(tissue), val(run_accession), path("*.bam")
    script:
    def starTmpDir =  "${params.outDir}/${taxon_id}/${run_accession}/alignment/tmp"
    def outFileNamePrefix = "${params.outDir}/${taxon_id}/${run_accession}/alignment/${run_accession}_"
    """
    rm -rf ${starTmpDir}
    STAR 
    --runThreadN ${task.cpus} \
    --twopassMode Basic 
    --runMode alignReads 
    --genomeDir ${file(genomeDir)} \
    --readFilesIn ${pair1} ${pair2} 
    --outFileNamePrefix ${outFileNamePrefix} \
    --readFilesCommand zcat 
    --outSAMattrRGline "ID:${run_accession}" 
    --outTmpDir ${starTmpDir} 
    --outSAMtype BAM SortedByCoordinate  

    --limitSjdbInsertNsj 2000000 \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --outSAMstrandField intronMotif 
    
    """
    
}
/*
STAR \
  --runThreadN 
  --genomeDir
  --readFilesIn
  --outSAMtype BAM SortedByCoordinate \
  --outSAMstrandField intronMotif \  It adds a tag (XS:A:+ or XS:A:-) tag to each alignment based on the strand of the splice junction. Used by StringTie
  --twopassMode Basic   Improves splice junction sensitivity/precision. Allows STAR to learn junctions in the first pass and improve mapping in the second.
  --outFilterIntronMotifs RemoveNoncanonicalUnannotated   Filters out introns that are non-canonical (non-GT/AG) and not annotated. 
  --limitSjdbInsertNsj 2000000 \  Limits the number of splice junctions to be inserted into the genome index. This is useful for large genomes or when there are many splice junctions.

to add??
  --outFilterType BySJout \  Use spliced junctions to filter alignments use known junctions-keep only those reads that contain junctions that passed filtering into SJ.out.tab
  --alignIntronMax 100000 \ filter long spurious introns

 
*/

  

