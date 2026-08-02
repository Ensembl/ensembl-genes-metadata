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
    * DOWNLOAD_FASTQS
    *
    * Download FASTQ files from ENA using the download_fastq.py script.
    *
    * Input:
    *   - meta: metadata map containing taxon_id, run_accession, url1, md5_1, url2, md5_2, paired, etc.
    *
    * Output:
    *   - FASTQ files (_1.fastq.gz and optionally _2.fastq.gz)
    *   - Software versions
    *
    * The module checks if the FASTQ files already exist before downloading to avoid redundant downloads.
    */

process DOWNLOAD_FASTQS {
    label "python"
    tag "${meta.taxonId}:${meta.run_accession}"
    maxForks 25
    //storeDir "${params.outDir}/$taxon_id/$run_accession"
    afterScript "sleep $params.files_latency"  // Needed because of file system latency
    //conda "$projectDir/pipelines/nextflow/modules/download_fastqs/environment.yml"

    input:
    tuple val(meta)

    //tuple val(taxon_id), val(gca), val(platform), val(paired), val(tissue), val(run_accession), val(genomeDir),  val(url1), val(md5_1), val(url2),  val(md5_2)

    output:
    tuple val(meta), path("*_1.fastq.gz"), path("*_2.fastq.gz", optional: true) , emit: fastq_file_output
    path "versions.yml", emit: versions_file


    script: 
    
    def fastq1 = "${params.outDir}/${meta.taxon_id}/${meta.run_accession}/${meta.run_accession}_1.fastq.gz"
    def fastq2 = meta.paired ? "${params.outDir}/${meta.taxon_id}/${meta.run_accession}/${meta.run_accession}_2.fastq.gz" : null

    def optionalArgs = meta.paired ? "--url2 ${meta.url2} --md5_2 ${meta.md5_2} --paired" : ""
    """
    if [ ! -s "${fastq1}" ] || { ${meta.paired} && [ ! -s "${fastq2}" ]; }; then 
    download_fastq.py \
        --taxon_id ${meta.taxon_id} \
        --run_accession ${meta.run_accession} \
        --url1 ${meta.url1} \
        --md5_1 ${meta.md5_1} \
        ${optionalArgs} \
        --outDir ${params.outDir} 

        ln -s ${params.outDir}/${meta.taxon_id}/${meta.run_accession}/*_1.fastq.gz ./
        ${meta.paired ? "ln -s ${params.outDir}/${meta.taxon_id}/${meta.run_accession}/*_2.fastq.gz ./" : ""}
        else
        echo "skipping"
        ln -s ${params.outDir}/${meta.taxon_id}/${meta.run_accession}/*_1.fastq.gz ./
        ${meta.paired ? "ln -s ${params.outDir}/${meta.taxon_id}/${meta.run_accession}/*_2.fastq.gz ./" : ""}
        fi

        
        cat <<-END_VERSIONS > versions.yml
            "${task.process}":
            download_fastq.py: \$(download_fastq.py --version 2>&1)
            python: \$(python --version | sed 's/Python //')
        END_VERSIONS
        """    
        }
