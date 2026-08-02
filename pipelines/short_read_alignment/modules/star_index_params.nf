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
    * STAR_INDEX_PARAMS
    *
    * Generate genome statistics using the star_genome_stats.py script.
    *
    * Input:
    *   - meta: metadata map containing taxon_id, gca, fasta_file, etc.
    *
    * Output:
    *   - Genome statistics in JSON format (stats.json)
    *   - Software versions
    *
    * The module uses the star_genome_stats.py script to generate genome statistics
    * and creates a versions.yml file with the software versions used.
    */
process STAR_INDEX_PARAMS {
    label 'python'
    tag "${meta.taxon_id}:${meta.gca}"

    input:
    val(meta)
    
    //when: meta.platform?.toString()?.toLowerCase() == 'illumina'

    output:
    //tuple val(taxon_id), val(genomeDir), val(platform),  val(tissue), val(run_accession), val(pair1), val(pair2)
    tuple val(meta), path("stats.json"), emit: genome_stats_output
    path "versions.yml", emit: versions_file

    script:
    """
    genome_stats.py ${meta.fasta_file} > stats.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star_genome_stats.py: $(star_genome_stats.py --version | awk '{print $2}')
        python: §$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}


