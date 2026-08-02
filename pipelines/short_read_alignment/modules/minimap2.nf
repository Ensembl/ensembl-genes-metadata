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
 * Align long-read sequencing data to a reference genome using Minimap2.
 *
 * The process maps Oxford Nanopore (ONT) and PacBio reads against a
 * pre-built Minimap2 genome index (.mmi). The alignment profile is
 * selected automatically based on the sequencing platform specified
 * in the metadata.
 */
process MINIMAP2 {
    tag "$meta.run_accession"
    label 'minimap2'
    storeDir "${params.outDir}/$meta.taxon_id/$meta.run_accession/alignment/"
    afterScript "sleep $params.files_latency"  // Needed because of file system latency

    input:
    //tuple val(taxon_id), val(genomeDir), val(platform), val(tissue), val(run_accession), val(input_file), path(minimap_index_file)
    tuple val(meta), path("${params.outDir}/$meta.taxon_id/$meta.run_accession/alignment/") ,path(minimap_index_file)

    output:
    //tuple val(taxon_id), val(genomeDir), val(tissue),val(platform),  val(run_accession), path("*.sam")
    tuple val(meta), path("*.sam") , emit: minimap_alignment
    path "versions.yml", emit: versions_file

    script:
    def sam_file = "${meta.run_accession}.sam"
    def profile = 
        "${meta.platform}" == 'ONT' ? '-x splice' :
        "${meta.platform}" == 'PacBio' ? '-x splice:hq' :
        '-ax splice'  // fallback
    //--secondary=no
    """
    minimap2 ${profile} -a -G ${params.max_intron_size} \
    --cs -N 1 -t ${params.cpus}  -u b ${minimap_index_file} \
    ${meta.pair1_path} -o ${sam_file}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
    END_VERSIONS
    """
}


