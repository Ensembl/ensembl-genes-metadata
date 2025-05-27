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
    storeDir "${params.outDir}/$taxon_id/$outputDir/alignment"
    afterScript "sleep $params.files_latency"  // Needed because of file system latency


    input:
    //tuple val(taxon_id), val(genomeDir), val(gca), val(platform), val(paired), val(tissue), val(run_accession), path(aligned_file)
    tuple val(taxon_id), val(genomeDir), val(tissue), val(outputDir),path(bamFiles)
    output:
    tuple val(taxon_id), val(genomeDir), val(tissue), val(outputDir),path("${tissue}.bam")

    script:
    """
    samtools merge ${tissue}.bam ${bamFiles.join(' ')}
    samtools index ${tissue}.bam ${tissue}.bam.bai
    """
}



