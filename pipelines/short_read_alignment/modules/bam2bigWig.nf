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
    * BAM2BIGWIG
    *
    * Convert a BAM file to BigWig format using bamCoverage.
    *
    * Input:
    *   - meta: metadata map containing taxon_id, platform, tissue, etc.
    *   - bam_file1: BAM file to convert (forward strand)
    *   - bam_file2: BAM file to convert (reverse strand, optional)
    *
    * Output:
    *   - BigWig files for forward and reverse strands
    *
    * The process uses bamCoverage to generate BigWig files from the input BAM files.
    */
process BAM2BIGWIG {
    tag "$bam_file1"
    label 'bamCoverage'
    publishDir "${meta.output_dir}", mode: "copy"
    afterScript "sleep $params.files_latency"  // Needed because of file system latency

    input:
    //tuple val(taxon_id), val(genomeDir), val(tissue), val(platform), val(output_dir), path(bam_file1),  path(bam_file2)
    tuple val(meta), path(bam_file1),  path(bam_file2)

    output:
    path("versions.yml"), emit: versions_file
    
    script:
    def bam_basename = bam_file1.baseName  // strips .bam 
    //def bam2_provided = bam_file2 ? true : false
    """
    ln -s ${meta.output_dir}/${bam_file1}.* .
    bamCoverage -b ${meta.output_dir}/${bam_file1} -o ${meta.output_dir}/${bam_basename}.bw --binSize 1 --numberOfProcessors ${task.cpus} 
    ln -s ${meta.output_dir}/${bam_basename}.bw .
    if [  -s "${bam_file2}" ]; then
      ln -s ${meta.output_dir}/${bam_file2}.* .
      bamCoverage -b ${meta.output_dir}/${bam_file2} -o ${meta.output_dir}/${bam_file2.baseName}.bw --binSize 1 --numberOfProcessors ${task.cpus}
      ln -s ${meta.output_dir}/${bam_file2.baseName}.bw .
    else
      echo "Reverse strand BAM not found, skipping..."
    fi 
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(bamCoverage --version 2>&1 | sed 's/^bamCoverage //')
    END_VERSIONS
    """
}


