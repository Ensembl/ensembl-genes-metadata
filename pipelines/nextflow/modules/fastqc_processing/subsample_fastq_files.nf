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
    label 'python'
    tag "subsampling $run_accession"
    storeDir "${params.outDir}/$taxon_id/$run_accession"
    afterScript "sleep $params.files_latency"  // Needed because of file system latency
    cache false 
    input:
    tuple val(taxon_id), val(gca), val(run_accession)
    
    //tuple val(pair1), val(pair2)

    output:
    //val(subsample_OUT)
    tuple val(taxon_id), val(gca), val(run_accession), path("*_1.fastq.gz.sub"),path("*_2.fastq.gz.sub")
    //val("${params.outDir}/${taxon_id}/${run_accession}/${run_accession}_1.fastq.gz.sub"), \
    //val("${params.outDir}/${taxon_id}/${run_accession}/${run_accession}_2.fastq.gz.sub")
    //subsampled_fastq_files = [Path(f"{fastq_file_1}.sub"), Path(f"{fastq_file_2}.sub")]

    script:
    def enscode1="/nfs/production/flicek/ensembl/genebuild/ftricomi/test_mypy/"
    def output_dir = "${params.outDir}/${taxon_id}/${run_accession}/"
    def fastq1 = "${output_dir}${run_accession}_1.fastq.gz"
    def fastq2 = "${output_dir}${run_accession}_2.fastq.gz"
    def subsampled_fastq1 = "${output_dir}${run_accession}_1.fastq.gz.sub"
    def subsampled_fastq2 = "${output_dir}${run_accession}_2.fastq.gz.sub"
    subsample_OUT=[]
    subsample_OUT.add([taxon_id:taxon_id, gca:gca, run_accession:run_accession, pair1:["${params.outDir}",taxon_id,run_accession,"${run_accession}_1.fastq.gz.sub"].join("/"), pair2:["${params.outDir}",taxon_id,run_accession,"${run_accession}_2.fastq.gz.sub"].join("/")])

    """
    chmod +x ${enscode1}/ensembl-anno/src/python/ensembl/tools/anno/transcriptomic_annotation/star.py
    if [ ! -s ${subsampled_fastq1} ] && [ ! -s ${subsampled_fastq2} ]; then 
       python3 ${enscode1}/ensembl-anno/src/python/ensembl/tools/anno/transcriptomic_annotation/star.py \
        --run_subsampling True --paired_file_1 ${fastq1} --paired_file_2 ${fastq2} --sampling_via_read_limit_percentage True \
        --subsample_read_limit 100000 --subsample_percentage 0.05 --num_threads 2  --output_dir ${output_dir};
    fi
    cp  ${output_dir}*.sub .
    echo '${subsample_OUT.toString()}'
    """

    //subsample_OUT.add([taxon_id:taxon_id, gca:gca, run_accession:run_accession, pair1:["${params.outDir}",taxon_id,run_accession,"${run_accession}_1.fastq.gz.sub"].join("/"), pair2:["${params.outDir}",taxon_id,run_accession,"${run_accession}_2.fastq.gz.sub"].join("/")])
}
    //"""
    //cp  ${output_dir}*.sub .
    //echo '${subsample_OUT.toString()}'
    //"""
