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
 * Generate a Minimap2 genome index (.mmi) from a reference FASTA.
 *
 * The process creates a Minimap2 index for long-read alignment.
 * If a valid .mmi index already exists in the genome directory,
 * the indexing step is skipped.
 */
process MINIMAP2_INDEX_GENOME {
    label 'minimap2'
    tag "${meta.taxon_id}:${meta.gca}"
    publishDir "${meta.genome_dir}", mode: 'copy'
    afterScript "sleep ${params.files_latency}"
    // Needed because of file system latency
    maxForks 10

    input:
    //tuple val(taxon_id), val(genomeDir), val(gca), val(platform), val(paired), val(tissue), val(run_accession), val(pair1)
    val(meta)

    output:
    tuple val(meta), path("*mmi"), emit: minimap_index
    path "versions.yml", emit: versions_file


    script:
    //def genomeDirPath= new File(genomeDir)
    //def fnaFiles = genomeDirPath.listFiles()?.findAll { it.name.endsWith('.fna') }
    if (!meta.fasta_file || meta.fasta_file.size() != 1) {
        throw new IllegalStateException("Expected exactly one .fna file in the directory: ${meta.fasta_file.parent}")
    }
    //def genomefilePath = fnaFiles[0]
    """
    if [ -z "\$(find "${meta.genome_dir}" -name '*.mmi' -type f -size +0c)" ]; then
        minimap2 --threads ${task.cpus} \
            -d ${meta.genome_dir}/genome.mmi ${meta.fasta_file}
    else
        echo "Minimap indexed genome already exists, skipping Minimap2 genomeGenerate step."
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
    END_VERSIONS

    """
}
