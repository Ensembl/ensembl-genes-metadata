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


process INDEX_GENOME {
    label 'star'
    tag "$run_accession"
    publishDir "${params.outDir}/$taxon_id/$gca/ncbi_dataset/", mode: 'copy'
    afterScript "sleep $params.files_latency"  // Needed because of file system latency

    input:
    tuple val(taxon_id), val(gca), val(run_accession), val(par_1), val(par_2), val(genome_file)

    output:
    tuple val(taxon_id), val(gca), val(run_accession), val(par_1), val(par_2), val("${params.outDir}/$taxon_id/$gca/ncbi_dataset/")
//    val("${params.outDir}/${taxon_id}/${run_accession}/${run_accession}_1.fastq.gz.sub"),\
//    val("${params.outDir}/${taxon_id}/${run_accession}/${run_accession}_2.fastq.gz.sub"),
//    val("${params.outDir}/$taxon_id/$gca/ncbi_dataset/")
    //subsampled_fastq_files = [Path(f"{fastq_file_1}.sub"), Path(f"{fastq_file_2}.sub")]

    script:

    def genomeDir = "${params.outDir}/$taxon_id/$gca/ncbi_dataset/"
    """
    if [ ! -s "${params.outDir}/${taxon_id}/${gca}/ncbi_dataset/Genome" ]; \
    then STAR --runThreadN ${task.cpus} --runMode genomeGenerate \
    --outFileNamePrefix ${genomeDir} --genomeDir ${genomeDir} \
    --genomeFastaFiles  ${file(genome_file)};fi
    
    """
}
/*
STAR --runMode genomeGenerate --genomeDir ./${assembly}_star_index --genomeFastaFiles ${genome} --runThreadN 4
if [ ! -e "/genome_dumps/Genome" ]; then /STAR --runThreadN 12 --runMode genomeGenerate 
--outFileNamePrefix /GCA_018555375.3//genome_dumps --genomeDir GCA_018555375.3//genome_dumps 
--genomeFastaFiles /GCA_018555375.3//genome_dumps/anguilla_rostrata_toplevel.fa;fi
*/
