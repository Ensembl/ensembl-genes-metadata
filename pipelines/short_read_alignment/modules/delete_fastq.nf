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
 * Delete downloaded FASTQ files after a successful alignment.
 *
 * This process removes compressed FASTQ (.fastq.gz) files from the sample
 * directory to free disk space. It assumes the alignment has completed
 * successfully and the BAM file has been generated.
 */
process DELETE_FASTQ {
    tag "$meta.run_accession"
    label 'default'

    input:
    //tuple val(taxon_id), val(genomeDir),  val(tissue),  path(aligned_file)
    //tuple val(taxon_id), val(genomeDir), val(tissue),val(platform),  val(output_dir), path(aligned_file)
    //tuple val(taxon_id), val(genomeDir), val(tissue), val(platform), val(run_accession), path(aligned_file)
    tuple val(meta), path(aligned_file)
    output:
    //tuple val(taxon_id), val(genomeDir), val(tissue), val(platform), val(run_accession), path(aligned_file)
    //tuple val(taxon_id), val(genomeDir), val(tissue), val(platform), val(output_dir), path(aligned_file)
    tuple val(meta), path(aligned_file) , emit:aligned_output
    path "versions.yml", emit: versions_file


    script:
    """
    rm -f ${params.outDir}/${meta.taxon_id}/${meta.run_accession}/*.gz
    echo "${meta.run_accession} deleted"
        cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: \$(rm --version | head -n1 | awk '{print \$NF}')
    END_VERSIONS
    """
    
}
