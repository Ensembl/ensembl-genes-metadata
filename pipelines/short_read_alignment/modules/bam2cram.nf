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
    * BAM2CRAM
    *
    * Convert aligned BAM files to CRAM format using samtools.
    *
    * Input:
    *   - meta: metadata map containing taxon_id, genomeDir, tissue, platform, output_dir, etc.
    *   - aligned_file: path to the input BAM file
    *
    * Output:
    *   - CRAM file (.cram)
    *   - Software versions
    *
    * The module uses samtools view to convert BAM to CRAM format and creates a symbolic link to the output CRAM file.
    */
process BAM2CRAM {
    tag "$aligned_file"
    label 'samtools'
    publishDir "${meta.output_dir}", mode: "copy"
    afterScript "sleep $params.files_latency"  // Needed because of file system latency

    input:
    //tuple val(taxon_id), val(genomeDir), val(tissue), val(platform), val(output_dir), path(aligned_file)
    tuple val(meta), path(aligned_file)
    output:
    tuple val(meta), path("*.cram"), emit:cram_output
    val("versions.yml"), emit: versions_file

    script:

    //def genomeDirPath= new File(genomeDir)
    //def genomeFile = genomeDirPath.listFiles().find { it.name.endsWith('fna') }
    def bam_basename = aligned_file.baseName  // strips .bam
    """
    samtools view -@ ${task.cpus}  -C -T ${meta.fasta_file} -o ${meta.output_dir}/${bam_basename}.cram ${aligned_file}
    ln -s ${meta.output_dir}/${bam_basename}.cram .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | awk '{print \$2}')
    END_VERSIONS
    """
    
}


