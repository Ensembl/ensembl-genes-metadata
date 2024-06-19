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
import java.nio.file.*

process INDEX_GENOME {
    label 'star'
    tag "$taxon_id:$run_accession:$gca"
    publishDir "${params.outDir}/$taxon_id/$gca/ncbi_dataset/", mode: 'copy'
    afterScript "sleep $params.files_latency"  // Needed because of file system latency
    maxForks 1
    input:
    tuple val(taxon_id), val(gca), val(run_accession), val(pair1), val(pair2), val(genomeDir)
    output:
    tuple val(taxon_id), val(gca), val(run_accession), val(pair1), val(pair2), val(genomeDir)
    //    val("${params.outDir}/${taxon_id}/${run_accession}/${run_accession}_1.fastq.gz.sub"),\
//    val("${params.outDir}/${taxon_id}/${run_accession}/${run_accession}_2.fastq.gz.sub"),
//    val("${params.outDir}/$taxon_id/$gca/ncbi_dataset/")
    //subsampled_fastq_files = [Path(f"{fastq_file_1}.sub"), Path(f"{fastq_file_2}.sub")]

    

    script:
    def d= new File("${genomeDir}")
    def genomefil=d.listFiles().find { it.name.endsWith('.fna') }
    def genomeFile=genomefil.absolutePath
    //genomeDirCopy=genomeDir
    //genomeFile = genomeDirCopy.map{file -> 
    //genomeFile = Channel.fromPath("${genomeDir}/*.fna")
   // def genomeDirPath = file(genomeDir)
    //def genomeFile = genomeDirPath.listFiles().find { it.name.endsWith('.fna') }
    //def genomeDir = "${params.outDir}/${taxon_id}/${gca}/ncbi_dataset/"
    // Find the first .fna file and assign it to genomeFile variable
    //def genomeFile = """find ${genomeDir} -maxdepth 1 -type f -name '*.fna' | head -n 1""".execute().text.trim()
    
    //genomeFile=\$( find ${genomeDir} -type f -name '*.fna')
    """
    if [ ! -s "${genomeDir}/Genome" ]; \
    then  rm -rf ${genomeDir}/_STARtmp ; STAR --runThreadN ${task.cpus} --runMode genomeGenerate \
    --outFileNamePrefix ${genomeDir} --genomeDir ${genomeDir} \
    --genomeFastaFiles  $genomeFile --outTmpDir _STARtmp;fi
    
    """
}
/*
STAR --runMode genomeGenerate --genomeDir ./${assembly}_star_index --genomeFastaFiles ${genome} --runThreadN 4
if [ ! -e "/genome_dumps/Genome" ]; then /STAR --runThreadN 12 --runMode genomeGenerate 
--outFileNamePrefix /GCA_018555375.3//genome_dumps --genomeDir GCA_018555375.3//genome_dumps 
--genomeFastaFiles /GCA_018555375.3//genome_dumps/anguilla_rostrata_toplevel.fa;fi
*/
