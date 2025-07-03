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

process INDEXING_FILES {
    tag "$aligned_file"
    label 'samtools'
    publishDir "${params.outDir}/$taxon_id/$output_dir/alignment", mode: 'copy'
    afterScript "sleep $params.files_latency"  // Needed because of file system latency
//
    input:
    //tuple val(taxon_id), val(genomeDir),  val(tissue),  path(aligned_file)
    tuple val(taxon_id), val(genomeDir), val(tissue),val(platform),  val(output_dir), path(aligned_file)
    val estension
    output:
    tuple val(taxon_id), val(genomeDir), val(tissue), val(platform), val(output_dir), path(aligned_file)



    script:
    def output_dir="${params.outDir}/$taxon_id/$output_dir/alignment"
    """
    if [ ! -s "${output_dir}/${aligned_file}.${estension}" ]; then
    samtools index ${output_dir}/${aligned_file} ${output_dir}/${aligned_file}.${estension}
    echo "${output_dir}/${aligned_file}.${estension}"
    else
    echo "skip file exists"
    fi
    """
    
}
