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


process SUBSAMPLE_FASTQ_FILES {
    label 'seqtk'
    tag "$taxon_id:$run_accession"
    storeDir "${params.outDir}/$taxon_id/$run_accession"
    afterScript "sleep $params.files_latency"  // Needed because of file system latency

    input:
    tuple val(taxon_id), val(gca), val(run_accession)

    output:
    tuple val(taxon_id), val(gca), val(run_accession), path("*_1.fastq.gz.sub"),path("*_2.fastq.gz.sub")

    script:
    def output_dir = "${params.outDir}/${taxon_id}/${run_accession}/"
    def fastq1 = "${output_dir}${run_accession}_1.fastq.gz"
    def fastq2 = "${output_dir}${run_accession}_2.fastq.gz"
    def subsampled_fastq1 = "${output_dir}${run_accession}_1.fastq.gz.sub"
    def subsampled_fastq2 = "${output_dir}${run_accession}_2.fastq.gz.sub"
    subsample_OUT=[]
    subsample_OUT.add([taxon_id:taxon_id, gca:gca, run_accession:run_accession, pair1:["${params.outDir}",taxon_id,run_accession,"${run_accession}_1.fastq.gz.sub"].join("/"), pair2:["${params.outDir}",taxon_id,run_accession,"${run_accession}_2.fastq.gz.sub"].join("/")])

    """
    seqtk sample -s100 ${fastq1} 50000 > ${subsampled_fastq1}
    seqtk sample -s100 ${fastq2} 50000 > ${subsampled_fastq2}
    cp  ${output_dir}*.sub .
    rm ${fastq1}
    rm ${fastq2}
    echo '${subsample_OUT.toString()}'
    """
}

