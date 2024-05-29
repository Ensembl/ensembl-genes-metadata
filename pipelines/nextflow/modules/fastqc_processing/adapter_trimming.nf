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

// module description 
process ADAPTER_TRIMMING {
    label 'python'
    tag "trimming $run_accession"
    publishDir "${params.outDir}/$taxon_id/$run_accession", mode: 'copy'
    input:
    tuple val(taxon_id), val(gca), val(run_accession),\
    val("${params.outDir}/${taxon_id}/${run_accession}/${run_accession}_1.fastq.gz.sub"), \
    val("${params.outDir}/${taxon_id}/${run_accession}/${run_accession}_2.fastq.gz.sub"), \
    val(overrepresented_sequences)
    //subsampled_fastq_files = [Path(f"{fastq_file_1}.sub"), Path(f"{fastq_file_2}.sub")]


    output:
    tuple val(taxon_id), val(gca), val(run_accession),\
    val("${params.outDir}/${taxon_id}/${run_accession}/trim_galore_output/${run_accession}_1.fastq.gz.fq"),\
    val("${params.outDir}/${taxon_id}/${run_accession}/trim_galore_output/${run_accession}_2.fastq.gz.fq")
    //no need to emit the path because the subsampled files will be in the run_accession dir _1_1

    script:
    """
    chmod +x ${params.enscode}/ensembl-anno/src/python/ensembl/tools/anno/transcriptomic_annotation/star.py 
            
    python ${params.enscode}/ensembl-anno/src/python/ensembl/tools/anno/transcriptomic_annotation/star.py \
        --run_trimming True --output_dir ${params.outDir}/${taxon_id}/${run_accession} \
        --short_read_fastq_dir ${params.outDir}/${taxon_id}/${run_accession} \
        --delete_pre_trim_fastq True --trim_galore_bin trim_galore --num_threads 2    
    """
}
