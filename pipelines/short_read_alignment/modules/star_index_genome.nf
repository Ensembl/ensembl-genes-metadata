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
    * STAR_INDEX_GENOME
    *
    * Generate STAR genome index for RNA-seq alignment.
    *
    * Input:
    *   - meta: metadata map containing taxon_id, gca, fasta_file, etc.
    *   - statsJson: path to the JSON file containing genome statistics
    *
    * Output:
    *   - Genome index directory
    *   - Software versions
    *
    * The module uses STAR to generate the genome index and creates a symbolic link to the output index directory.
    */  
process STAR_INDEX_GENOME {
    label 'star'
    tag "${meta.taxon_id}:${meta.gca}"
    publishDir "${meta.genome_dir}", mode: 'copy'
    afterScript "sleep $params.files_latency"  // Needed because of file system latency
    maxForks 1

    input:
    //tuple val(taxon_id), val(genomeDir), val(gca), val(platform), val(paired), val(tissue), val(run_accession), val(pair1), val(pair2)
    tuple val(meta),path(statsJson)



    output:
    //tuple val(taxon_id), val(genomeDir), val(platform),  val(tissue), val(run_accession), val(pair1), val(pair2)
    val(meta), emit: genome_index_output
    path "versions.yml", emit: versions_file
    
    //when: meta.platform?.toString()?.toLowerCase() == 'illumina'

    script:
    //def stats = new groovy.json.JsonSlurper().parse(file('stats.json'))
    def genomeDir = meta.genome_dir
    def limitBAMsortRAM = (task.memory.toBytes() * 0.8) as long
    //def genomeDirPath= new File(genomeDir)
    //def genomeIndexFile = genomeDirPath.listFiles()?.find { it.name.endsWith('Genome') }
    //log.info("Genome index file: ${genomeIndexFile?.absolutePath}")
    
    //def genomefilePath = genomeDirPath.listFiles()?.find { it.name.endsWith('.fna') }

    
    """
    genomeSAindexNbases=\$(sed -n 's/.*"genomeSAindexNbases"[[:space:]]*:[[:space:]]*\\([0-9][0-9]*\\).*/\\1/p' "${statsJson}")
    genomeChrBinNbits=\$(sed -n 's/.*"genomeChrBinNbits"[[:space:]]*:[[:space:]]*\\([0-9][0-9]*\\).*/\\1/p' "${statsJson}")
    echo \$genomeSAindexNbases
    echo \$genomeChrBinNbits
    if [ ! -s "${genomeDir}/Genome" ]; then
        rm -rf ${genomeDir}/_STARtmp && \
        STAR --runThreadN ${task.cpus} --runMode genomeGenerate \
            --outFileNamePrefix ${genomeDir} \
            --genomeDir ${genomeDir} \
            --genomeSAindexNbases "\$genomeSAindexNbases" \
            --genomeChrBinNbits "\$genomeChrBinNbits" \
            --genomeFastaFiles "${meta.fasta_file}" \
            --outTmpDir _STARtmp \
            --limitBAMsortRAM ${limitBAMsortRAM};   
    
    
    else 
    
    echo "Genome index already exists, skipping STAR genomeGenerate step."
    
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    star : \$(STAR --version 2>&1)
    END_VERSIONS  
    """ 
}


