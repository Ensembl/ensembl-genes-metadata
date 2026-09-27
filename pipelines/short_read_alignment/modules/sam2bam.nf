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
 * Convert a SAM alignment file to BAM format and generate a BAM index.
 *
 * The process uses samtools to convert the SAM output produced by
 * Minimap2 into a compressed BAM file and creates the corresponding
 * BAM index (.bai) for downstream analyses.
 */

process SAM2BAM {
    tag "${meta.run_accession}"
    label 'samtools'
    storeDir "${meta.alignment_dir}"
    afterScript "sleep $params.files_latency"  // Needed because of file system latency

    input:
    tuple val(meta), path(sam_file)


    output:
    tuple val(meta), path("*.bam"), emit:sam_output
    path "versions.yml", emit: versions_file

    script:
    """
    samtools sort \
    -@ ${task.cpus} \
    -o ${meta.run_accession}.bam \
    ${sam_file}

    samtools index \
    -@ ${task.cpus} \
    ${meta.run_accession}.bam
    
    SAMTOOLS_VERSION=\$(samtools --version | head -n1)
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \${SAMTOOLS_VERSION}
    END_VERSIONS
    """
}


