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
The process
*/
process INDEXING_FILES {
    tag "${meta.run_accession}"
    label 'samtools'
    publishDir "${params.outDir}/${meta.taxon_id}/${meta.output_dir}/alignment", mode: 'copy'
    afterScript "sleep $params.files_latency"  // Needed because of file system latency
//
    input:
    //tuple val(taxon_id), val(genomeDir),  val(tissue),  path(aligned_file)
    //tuple val(taxon_id), val(genomeDir), val(tissue),val(platform),  val(output_dir), path(aligned_file)
    tuple val(meta), path(aligned_file)
    val extension
    output:
    //tuple val(taxon_id), val(genomeDir), val(tissue), val(platform), val(output_dir), path(aligned_file)
    tuple val(meta), path(aligned_file), emit:aligned_output
    path "versions.yml", emit: versions_file


    script:
    //def output_dir="${params.outDir}/${meta.taxon_id}/${meta.output_dir}/alignment"
    """
    if [ ! -s "${meta.output_dir}/${aligned_file}.${extension}" ] || [ ! -s "${meta.output_dir}/${aligned_file}.csi" ]; then
    samtools index -c  \
    ${meta.output_dir}/${aligned_file} ${meta.output_dir}/${aligned_file}.${extension} \
    -@ ${task.cpus}
    echo "${meta.output_dir}/${aligned_file}.${extension}"
    else
    echo "skip file exists"
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | awk '{print \$2}')
    END_VERSIONS
    """
    
}
