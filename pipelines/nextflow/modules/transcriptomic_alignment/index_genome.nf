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
    tag "$taxon_id:$gca"
    publishDir "${params.outDir}/$taxon_id/$gca/ncbi_dataset/", mode: 'copy'
    afterScript "sleep $params.files_latency"  // Needed because of file system latency
    maxForks 10

    input:
    tuple val(taxon_id), val(gca), val(run_accession), val(pair1), val(pair2), val(genomeDir)

    output:
    tuple val(taxon_id), val(gca), val(run_accession), val(pair1), val(pair2), val(genomeDir)

    script:
    def genomeDirPath= new File(genomeDir)
    def genomeIndexFile = genomeDirPath.listFiles().find { it.name.endsWith('Genome') }
    log.info("${genomeIndexFile}")
    if (!genomeIndexFile || genomeIndexFile.length() == 0) {
    //new File("${genomeDir}/Genome")
    //if (!genomeIndexFile.exists() || genomeIndexFile.length() == 0) {
    // Read the .fna file and perfor
    // d= new File("${genomeDir}")
    def genomefilePath=d.listFiles().find { it.name.endsWith('.fna') }
    def genomeFile=genomefilePath.absolutePath
    // Function to calculate the min value
    def min = { a, b -> a < b ? a : b }

    // Initialize variables
    def numberOfReferences = 0
    def genomeLength = 0
    def retryCount = 0
    def maxRetries = 3
    def filesValid = false
    while (!filesValid && retryCount < maxRetries) {
        // Read the FASTA file line by line
        genomefilePath.eachLine { line ->
            if (line.startsWith('>')) {
                numberOfReferences++
            } else {
            genomeLength += line.trim().length()
        }
    }
    log.info("numberOfReferences: ${numberOfReferences}")
    log.info("genomeLength: ${genomeLength.abs()}")
    // Read the FASTA file
    //def fastaContent = genomefilePath.text

    // Calculate the number of references
    //def numberOfReferences = fastaContent.count('>')

    // Calculate the genome length (excluding header lines)
    //def genomeLength = fastaContent.split('\n').findAll { !it.startsWith('>') }.join('').length()

     if (genomeLength != 0) {
       filesValid = true
        }
    }
    // Define read length (you may need to adjust this based on your data)
    def readLength = 100

    // Calculate genomeSAindexNbases
    def genomeSAindexNbases = min(14, (Math.log(genomeLength.abs() as Double) / Math.log(2) / 2 - 1) as int)

    // Calculate genomeChrBinNbits
    //def genomeChrBinNbits = min(18, (Math.log(Math.max(genomeLength as Double/ numberOfReferences as Double, readLength as Double)) / Math.log(2)) as int)
    def genomeChrBinNbits = min(18, (Math.log(Math.max((genomeLength.abs() / numberOfReferences) as Double, readLength as Double)) / Math.log(2)) as int)

    // Print the calculated values for debugging
    log.info("genomeSAindexNbases: ${genomeSAindexNbases}")
    log.info("genomeChrBinNbits: ${genomeChrBinNbits}")

    // Execute STAR command with calculated parameters #if [ ! -s "${genomeDir}/Genome" ]; then \
    """
    if [ ! -s "${genomeDir}/Genome" ]; \
    then
    rm -rf ${genomeDir}/_STARtmp ; 
    STAR  --runThreadN ${task.cpus} --runMode genomeGenerate \
    --outFileNamePrefix ${genomeDir} --genomeDir ${genomeDir} \
    --genomeSAindexNbases ${genomeSAindexNbases} \
    --genomeChrBinNbits ${genomeChrBinNbits} \
    --genomeFastaFiles  ${genomeFile} --outTmpDir _STARtmp;fi
    """
    } else {
    """
    echo "Genome index already exists, skipping STAR genomeGenerate step."
    """
    }

}


