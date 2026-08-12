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

process WRITE_REPORT {
    label "python"
    tag "${meta.taxon_id}"
    //storeDir "${params.outDir}/$taxon_id/$run_accession"
    afterScript "sleep $params.files_latency"  // Needed because of file system latency

    input:
    val(meta)

    output:
    path "versions.yml", emit: versions_file


    script: 
    """
    write_report.py \
        --csv_path ${meta.csv_path} \
        --base_dir ${meta.output_dir} \
        --output_csv ${meta.output_dir}/report.csv \
        --merge_tissue ${params.mergeTissue} \
        --bam2bigWig ${params.bam2bigWig} \

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
        write_report.py: \$(write_report.py --version 2>&1)
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """    
}
