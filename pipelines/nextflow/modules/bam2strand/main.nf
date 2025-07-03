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

process BAM2STRAND {
    tag "$aligned_file"
    label 'samtools'
    storeDir "$output_dir"
    //publishDir "${params.outDir}/$taxon_id/$outut_dir/alignment", mode: 'copy'
    afterScript "sleep $params.files_latency"  // Needed because of file system latency

    input:
    tuple val(taxon_id), val(genomeDir), val(tissue), val(platform), val(output_dir), path(aligned_file)


    output:

    tuple val(taxon_id), val(genomeDir), val(tissue), val(platform), val(output_dir), path("*_forward_strand.bam"), path("*_reverse_strand.bam")



    script:
    def bam_basename = aligned_file.baseName  // strips .bam
    """
    # Plus strand
    samtools view -h ${aligned_file} | grep -E '^@|XS:A:\\+' | samtools view -Sb - > ${output_dir}/${bam_basename}_forward_strand.bam
    samtools index ${output_dir}/${bam_basename}_forward_strand.bam ${output_dir}/${bam_basename}_forward_strand.bam.bai

    # Minus strand
    samtools view -h ${aligned_file} | grep -E '^@|XS:A:-' | samtools view -Sb - > ${output_dir}/${bam_basename}_reverse_strand.bam
    samtools index ${output_dir}/${bam_basename}_reverse_strand.bam ${output_dir}/${bam_basename}_reverse_strand.bam.bai

    """
}


