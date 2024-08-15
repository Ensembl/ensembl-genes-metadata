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

process CONVERT_BAM_TO_CRAM {
    tag "$run_accession"
    label 'samtools'
    storeDir "${params.outDir}/$taxon_id/$run_accession/star/"
    afterScript "sleep $params.files_latency"  // Needed because of file system latency

    input:
    tuple val(taxon_id), val(gca), val(run_accession), val(genomeDir),val(bamFile)

    output:
    tuple val(taxon_id), val(gca), val(run_accession), path('*.cram')

    script:
    def genomeDirPath= new File(genomeDir)
    def genomeIndexFile = genomeDirPath.listFiles().find { it.name.endsWith('Genome') }
    """
    samtools view -C -T reference.fa -o output.cram input.bam

    """
    
}
