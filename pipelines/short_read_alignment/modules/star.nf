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

/* 
The STAR process is used to align short reads to a reference genome. 
STAR documentation https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf

STAR \
  --runThreadN 
  --genomeDir
  --readFilesIn
  --outSAMtype BAM SortedByCoordinate \
  --outSAMstrandField intronMotif \  It adds a tag (XS:A:+ or XS:A:-) tag to each alignment based on the strand of the splice junction. Used by StringTie
  --twopassMode Basic   Improves splice junction sensitivity/precision. Allows STAR to learn junctions in the first pass and improve mapping in the second.
  --outFilterIntronMotifs RemoveNoncanonicalUnannotated   Filters out introns that are non-canonical (non-GT/AG) and not annotated. 
  --limitSjdbInsertNsj 2000000 \  Limits the number of splice junctions to be inserted into the genome index. This is useful for large genomes or when there are many splice junctions.

Consider using the following parameters for STAR alignment:
  --outFilterType BySJout \  Use spliced junctions to filter alignments use known junctions-keep only those reads that contain junctions that passed filtering into SJ.out.tab
  --alignIntronMax 100000 \ filter long spurious introns

*/

process STAR {
    tag "$meta.run_accession"
    label 'star'
    publishDir "$meta.alignment_dir", mode: 'copy'
    afterScript "sleep $params.files_latency"  // Needed because of file system latency

    input:
    //tuple val(taxon_id), val(genomeDir), val(platform), val(tissue), val(run_accession), val(pair1), val(pair2)
    val(meta)

    output:
    //tuple val(taxon_id), val(genomeDir), val(gca), val(platform), val(paired), val(tissue), val(run_accession), path("*_Aligned.sortedByCoord.out.bam")
    //tuple val(taxon_id), val(genomeDir), val(tissue), val(platform), val(run_accession), path("*.bam")
    tuple val(meta), path("${meta.run_accession}_Aligned.sortedByCoord.out.bam"), emit: star_output
    path "versions.yml", emit: versions_file
    
    script:
    def starTmpDir =  "${meta.alignment_dir}/tmp"
    def outFileNamePrefix = "${meta.run_accession}_"
    //def outFileNamePrefix = "${params.outDir}/${meta.taxon_id}/${meta.run_accession}/alignmenti/${meta.run_accession}_"
    def limitBAMsortRAM = (task.memory.toBytes() * 0.9) as long
    """
    if [ ! -s "$meta.alignment_dir/${meta.run_accession}_Aligned.sortedByCoord.out.bam" ]; then
    rm -rf ${starTmpDir}
    mkdir -p ${meta.alignment_dir}
    STAR \
    --runThreadN ${task.cpus} \
    --twopassMode Basic \
    --runMode alignReads \
    --genomeDir ${meta.genome_dir} \
    --readFilesIn ${meta.fastq1} ${meta.fastq2} \
    --outFileNamePrefix ${outFileNamePrefix} \
    --readFilesCommand zcat \
    --outSAMattrRGline "ID:${meta.run_accession}" \
    --outTmpDir ./tmp \
    --outSAMtype BAM SortedByCoordinate  \
    --limitSjdbInsertNsj ${params.limitSjdbInsertNsj} \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --outSAMstrandField intronMotif \
    --limitBAMsortRAM ${limitBAMsortRAM} 
    
    else
    echo "skip file exists"
    ln -s ${meta.alignment_dir}/*.bam ./
    fi

    STAR_VERSION=\$(STAR --version | head -n1)
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \${STAR_VERSION}
    END_VERSIONS
    
    """
}    


