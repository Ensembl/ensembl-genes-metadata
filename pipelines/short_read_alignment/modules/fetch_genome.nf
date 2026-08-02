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
This process fetches the genome file for a given GCA accession. 
If a genome file is provided as part of the metadata, 
it uses that file instead of downloading it. 
The fetched genome file is saved with the name "genome.fna".
*/
process FETCH_GENOME {
    tag "${meta.gca}:genome"
    label 'fetch_file'
    storeDir "${params.outDir}/${meta.taxon_id}/${meta.gca}"
    afterScript "sleep $params.files_latency"  // Needed because of file system latency
    maxForks 1

    input:
    val meta
    //tuple val(taxon_id), val(gca), val(platform), val(paired), val(tissue), val(run_accession), val(url1), val(md5_1), val(url2),  val(md5_2) 
    
    output:
    tuple val(meta), path("*.fna"), emit: fasta_file_output
    path "versions.yml", emit: versions_file
    //tuple val(taxon_id), val(gca), val(platform), val(paired), val(tissue), val(run_accession), val("${params.outDir}/$taxon_id/$gca/"), val(url1), val(md5_1), val(url2), val(md5_2)
    
    script:
    """
    if [[ -f "${meta.genome_file}" ]]; then
            echo "Using provided genome file: ${meta.genome_file}"
            cp -L "${meta.genome_file}" genome.fna
        else
            fetch_genome.py \
                --output_dir . \
                --gca ${meta.gca} \
                --ncbi_base ${params.ncbiBaseUrl}

            downloaded_genome=\$(find . -maxdepth 1 -type f -name "*.fna" | head -n 1)

            if [[ -z "\$downloaded_genome" ]]; then
                echo "No genome FASTA found for ${meta.gca}" >&2
                exit 1
            fi
        fi
        
        # Create versions file

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fetch_genome.py: \$(fetch_genome.py --version 2>&1)
            python: \$(python --version | sed 's/Python //')
        END_VERSIONS
        
    """

    }
