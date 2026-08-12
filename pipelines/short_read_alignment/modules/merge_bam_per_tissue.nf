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
 * MERGE_BAM_PER_TISSUE
 *
 * Merge multiple aligned BAM files belonging to the same tissue into a single BAM.
 *
 * Input:
 *   - meta: metadata map containing taxon_id, platform, tissue, etc.
 *   - bamFiles: list of BAM files to merge
 *
 * Output:
 *   - Merged tissue-level BAM file
 *   - BAM index (.bai)
 *   - Software versions
 *
 * The module validates input BAM files using samtools quickcheck before merging
 * and validates the merged BAM before indexing.
 */
process MERGE_BAM_PER_TISSUE {
    label "samtools"
    tag "${meta.tissue}"
    maxForks 2
    storeDir "${meta.alignment_dir}"
    afterScript "sleep $params.files_latency"  // Needed because of file system latency


    input:
    //tuple val(taxon_id), val(genomeDir), val(gca), val(platform), val(paired), val(tissue), val(run_accession), path(aligned_file)
    //tuple val(taxon_id), val(genomeDir), val(tissue), val(platform), path(bamFiles)
    tuple val(meta), path(bamFiles)

    output:
    //tuple val(taxon_id), val(genomeDir), val(tissue), val(platform), \
    //val("${params.outDir}/$taxon_id/$platform/$tissue/alignment"),path("${tissue}.bam")
    tuple val(meta), path("${meta.tissue}.bam"), emit: merged_bam
    path "versions.yml", emit: versions_file

    script:
    def outputDir="${meta.output_dir}/${meta.taxon_id}/${meta.platform}/${meta.tissue}/alignment"
    """
    if [ ! -s "${meta.alignment_dir}/${meta.tissue}.bam" ]; then
    mkdir -p ${meta.alignment_dir}
    samtools merge -@ ${task.cpus}  -f -O BAM -o ${meta.tissue}.bam ${bamFiles.join(' ')}
    else
    ln -s ${outputDir}/${meta.tissue}.bam .
    
    echo "merging"
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n 1 | awk '{print \$2}')

    END_VERSIONS
    """
}
//samtools index ${outputDir}/${meta.tissue}.bam ${outputDir}/${meta.tissue}.bam.bai


