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

process MERGE_BAM_PER_TISSUE {
    label "samtools"
    tag "${taxon_id}"
    maxForks 25
    storeDir "${params.outDir}/$taxon_id/$platform/$tissue/alignment"
    afterScript "sleep $params.files_latency"  // Needed because of file system latency


    input:
    //tuple val(taxon_id), val(genomeDir), val(gca), val(platform), val(paired), val(tissue), val(run_accession), path(aligned_file)
    tuple val(taxon_id), val(genomeDir), val(tissue), val(platform), path(bamFiles)
    output:
    tuple val(taxon_id), val(genomeDir), val(tissue), val(platform), val("${params.outDir}/$taxon_id/$platform/$tissue/alignment"),path("${tissue}.bam")

    script:
    def outputDir="${params.outDir}/$taxon_id/$platform/$tissue/alignment"
    """
    mkdir -p ${outputDir}
    samtools merge -f -o ${outputDir}/${tissue}.bam ${bamFiles.join(' ')}
    samtools index ${outputDir}/${tissue}.bam ${outputDir}/${tissue}.bam.bai
    ln -s ${outputDir}/${tissue}.bam .
    """
}



