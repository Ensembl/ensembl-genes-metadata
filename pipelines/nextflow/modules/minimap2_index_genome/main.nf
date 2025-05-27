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


process MINIMAP2_INDEX_GENOME {
    label 'minimap2'
    tag "${taxon_id}:${gca}"
    publishDir "${genomeDir}", mode: 'copy'
    afterScript "sleep ${params.files_latency}"
    // Needed because of file system latency
    maxForks 10

    input:
    tuple val(taxon_id), val(genomeDir), val(gca), val(platform), val(paired), val(tissue), val(run_accession), val(pair1)


    output:
    tuple val(taxon_id), val(genomeDir), val(platform), val(tissue), val(run_accession), val(pair1), path("*mmi")


    script:
    def genomeDirPath= new File(genomeDir)
    def fnaFiles = genomeDirPath.listFiles()?.findAll { it.name.endsWith('.fna') }
    if (!fnaFiles || fnaFiles.size() != 1) {
        throw new IllegalStateException("Expected exactly one .fna file in the directory: ${genomeDirPath}, but found ${fnaFiles?.size() ?: 0}")
    }
    def genomefilePath = fnaFiles[0]
    """
    if [ ! -z "\$(find "${genomeDir}" -name '*.mmi' -type f -size +0c)" ]; then
        minimap2 --threads ${task.cpus} \
            -d ${genomefilePath}.mmi ${genomefilePath}
    else
        echo "Minimap indexed genome already exists, skipping Minimap2 genomeGenerate step."
    fi
    """
}