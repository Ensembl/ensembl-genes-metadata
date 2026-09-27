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
   * BAM2STRAND
   *
   * Split a BAM file into forward and reverse strand BAM files.
   *
   * Input:
   *   - meta: metadata map containing taxon_id, platform, tissue, etc.
   *   - aligned_file: BAM file to split
   *
   * Output:
   *   - Forward strand BAM file
   *   - Reverse strand BAM file
   *   - Software versions
   *
   * The process uses samtools to filter the input BAM file based on the strand information
   * and generates two separate BAM files for forward and reverse strands.
   */
process BAM2STRAND {
   tag "$aligned_file"
   label 'samtools'
   publishDir "${meta.alignment_dir}", mode: 'copy'
   afterScript "sleep $params.files_latency"  // Needed because of file system latency

   input:
   tuple val(meta), path(aligned_file)

   output:
   tuple val(meta), path("*_forward_strand.bam"), path("*_reverse_strand.bam"), path("*_forward_strand.bam.csi"), path("*_reverse_strand.bam.csi"), emit:aligned_output
   val(meta), emit:meta_value
   path "versions.yml", emit: versions_file

   script:
   def bam_basename = aligned_file.baseName  // strips .bam
   """
   # Plus strand
   # samtools view -h ${aligned_file} | grep -E '^@|XS:A:\\+' | samtools view -Sb - > ${bam_basename}_forward_strand.bam
   samtools view -b -f 0x2 -F 0x10 ${aligned_file} > ${bam_basename}_forward_strand.bam
   
   samtools index -c ${bam_basename}_forward_strand.bam ${bam_basename}_forward_strand.bam.csi

   # Minus strand
   #samtools view -h ${aligned_file} | grep -E '^@|XS:A:-' | samtools view -Sb - > ${bam_basename}_reverse_strand.bam
   samtools view -b -f 0x2 -f 0x10 ${aligned_file} > ${bam_basename}_reverse_strand.bam

   samtools index -c ${bam_basename}_reverse_strand.bam ${bam_basename}_reverse_strand.bam.csi

   # Optional sanity check
   if [ ! -s "${bam_basename}_forward_strand.bam" ]; then
      echo "Warning: forward strand BAM is empty"
      exit 1
   fi
   if [ ! -s "${bam_basename}_reverse_strand.bam" ]; then
      echo "Warning: reverse strand BAM is empty"
      exit 1
   fi

   SAMTOOLS_VERSION=\$(samtools --version | head -n1)
   cat <<-END_VERSIONS > versions.yml
   "${task.process}":
      samtools: \${SAMTOOLS_VERSION}
   END_VERSIONS
   """
}


