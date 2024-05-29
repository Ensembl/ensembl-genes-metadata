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
    publishDir "${params.outDir}/$taxon_id/$run_accession", mode: 'copy'
    afterScript "sleep $params.files_latency"  // Needed because of file system latency
    
    input:
    tuple val(taxon_id), val(gca), val(run_accession)
    
    tuple path(pair1), path(pair2)

    output:
    tuple val(taxon_id), val(gca), val(run_accession),\
    val("${params.outDir}/${taxon_id}/${run_accession}/${run_accession}_1.fastq.gz.sub"),\
    val("${params.outDir}/${taxon_id}/${run_accession}/${run_accession}_2.fastq.gz.sub")
    //subsampled_fastq_files = [Path(f"{fastq_file_1}.sub"), Path(f"{fastq_file_2}.sub")]

    script:
    """
    chmod +x ${params.enscode}/ensembl-anno/src/python/ensembl/tools/anno/transcriptomic_annotation/star.py 
            
    python ${params.enscode}/ensembl-anno/src/python/ensembl/tools/anno/transcriptomic_annotation/star.py \
        --run_subsampling True --paired_file_1 ${pair1} --paired_file_2 ${pair2} --sampling_via_read_limit_percentage\
        --subsample_read_limit 100000 --subsample_percentage 0.10 --num_threads 2   
    """
}
