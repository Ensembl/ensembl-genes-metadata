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


process STAR_INDEX_GENOME {
    label 'star'
    tag "$taxon_id:$gca"
    publishDir "${genomeDir}", mode: 'copy'
    afterScript "sleep $params.files_latency"  // Needed because of file system latency
    maxForks 10

    input:
    tuple val(taxon_id), val(genomeDir), val(gca), val(platform), val(paired), val(tissue), val(run_accession), val(pair1), val(pair2)


    output:
    tuple val(taxon_id), val(genomeDir), val(platform),  val(tissue), val(run_accession), val(pair1), val(pair2)


    script:
    def genomeDirPath= new File(genomeDir)
    def genomeIndexFile = genomeDirPath.listFiles()?.find { it.name.endsWith('Genome') }
    log.info("Genome index file: ${genomeIndexFile?.absolutePath}")
    if (!genomeIndexFile || genomeIndexFile.length() == 0) {
    //new File("${genomeDir}/Genome")
    //if (!genomeIndexFile.exists() || genomeIndexFile.length() == 0) {
    // Read the .fna file and perfor
    // d= new File("${genomeDir}")
    def genomefilePath = genomeDirPath.listFiles()?.find { it.name.endsWith('.fna') }
    if (!genomefilePath) {
        throw new IllegalStateException("No .fna file found in the directory: ${genomeDirPath}")
    }
    def genomeFile=genomefilePath.absolutePath
    // Function to calculate the min value
    def min = { a, b -> a < b ? a : b }

    // Initialize variables
    def numberOfReferences = 0
    def genomeLength = 0
    //def retryCount = 0
    def maxRetries = 3
    //def filesValid = false
    (0..<maxRetries).find { retryCount ->
    try {
    //while (!filesValid && retryCount < maxRetries) {
        // Read the FASTA file line by line
        genomefilePath.eachLine { line ->
            if (line.startsWith('>')) {
                    numberOfReferences += 1
            } else {
            genomeLength += line.trim().length()
        }
    }
    log.info("numberOfReferences: ${numberOfReferences}")
    log.info("genomeLength: ${genomeLength.abs()}")
    if (numberOfReferences == 0 || genomeLength == 0) {
                throw new IllegalStateException("The .fna file is empty or invalid: ${genomefilePath}")
            }
    return true        
    // Read the FASTA file
    //def fastaContent = genomefilePath.text

    // Calculate the number of references
    //def numberOfReferences = fastaContent.count('>')

    // Calculate the genome length (excluding header lines)
    //def genomeLength = fastaContent.split('\n').findAll { !it.startsWith('>') }.join('').length()

    //if (genomeLength > 0) {
    //    return true
    //}
    }
    catch (Exception e) {
        log.warn("Attempt ${retryCount + 1} failed: ${e.message}")
        return false
    }
    } ?: { throw new RuntimeException("File validation failed after ${maxRetries} retries.") }()
    // Define read length (you may need to adjust this based on your data)
    def readLength = 100

    // Calculate genomeSAindexNbases
    def genomeSAindexNbases = genomeLength > 0 ? 
        min(14, (Math.log(genomeLength.abs() as Double) / Math.log(2) / 2 - 1) as int) : 
        0 // Default value or handle the case where genomeLength is 0

    def genomeChrBinNbits = numberOfReferences > 0 ? 
        min(18, (Math.log(Math.max((genomeLength.abs() / numberOfReferences) as Double, readLength as Double)) / Math.log(2)) as int) : 
        0
    //def genomeChrBinNbits = min(18, (Math.log(Math.max(genomeLength as Double/ numberOfReferences as Double, readLength as Double)) / Math.log(2)) as int)
    //def genomeChrBinNbits = min(18, (Math.log(Math.max((genomeLength.abs() / numberOfReferences) as Double, readLength as Double)) / Math.log(2)) as int)

    // Print the calculated values for debugging
    log.info("genomeSAindexNbases: ${genomeSAindexNbases}")
    log.info("genomeChrBinNbits: ${genomeChrBinNbits}")
    """
    if [ ! -s "${genomeDir}/Genome" ]; then
        rm -rf ${genomeDir}/_STARtmp && \
        STAR --runThreadN ${task.cpus} --runMode genomeGenerate \
            --outFileNamePrefix ${genomeDir} \
            --genomeDir ${genomeDir} \
            --genomeSAindexNbases ${genomeSAindexNbases} \
            --genomeChrBinNbits ${genomeChrBinNbits} \
            --genomeFastaFiles ${genomeFile} \
            --outTmpDir _STARtmp
    fi
    """
    } else {
    """
    echo "Genome index already exists, skipping STAR genomeGenerate step."
    """
    }    
    }


