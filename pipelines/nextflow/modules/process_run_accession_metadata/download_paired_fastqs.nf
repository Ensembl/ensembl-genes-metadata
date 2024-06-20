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

include { getRunTable } from '../utils.nf'
include { getDataFromTable } from '../utils.nf'
import groovy.json.JsonSlurper
import org.apache.commons.codec.digest.DigestUtils
import java.nio.file.Files
import java.nio.file.Path

process DOWNLOAD_PAIRED_FASTQS {
    label "default"
    tag "${taxon_id}:${run_accession}"
    maxForks 10
    storeDir "${params.outDir}/$taxon_id/$run_accession"
    afterScript "sleep $params.files_latency"  // Needed because of file system latency
    
    input:
    tuple val(taxon_id), val(gca), val(run_accession)
    path(dataFileQuery)

    output:
    //val(fastqMetadata)
    val(download_OUT)
    //tuple (val(taxon_id), val(gca), val(run_accession), val("${params.outDir}/${taxon_id}/${run_accession}/${run_accession}_1.fastq.gz"), val("${params.outDir}/${taxon_id}/${run_accession}/${run_accession}_2.fastq.gz"),val("${params.outDir}/$taxon_id/$run_accession/${dataFileQuery}"))

    //when:
    //dataFiles.size() == 2 && qc_status != 'FILE_ISSUE'


    script:
    //def fastqMetadata = []
    //def fileContent = Files.readString(dataFileQuery)
    //def parsedJson = new JsonSlurper().parseText(fileContent)
    //def dataFiles = parsedJson.data_files
   //file(${params.outDir}/$taxon_id/$run_accession/dataFileQuery) 
    // Parse the JSON data
    def parsedJson = new JsonSlurper().parse(file("${params.outDir}/$taxon_id/$run_accession/$dataFileQuery"))
    def dataFiles = parsedJson.data_files

    // Check if there are exactly two data files
    if (dataFiles.size() != 2) {
        println "Expected two data files, but found ${dataFiles.size()}. Skipping the process."
        return []
    }

    def file1 = dataFiles[0]
    def file2 = dataFiles[1]
    println(file1)
    println(file2)
    // Extract URLs and MD5 checksums
    def url1 = file1.file_url
    def md5_1 = file1.md5
    def url2 = file2.file_url
    def md5_2 = file2.md5
    //def qc_status = getRunTable(run_accession, "qc_status")
    def qc_status = getDataFromTable("qc_status", "run", "run_accession", run_accession)[0].qc_status
    // Check for file issues and QC status
    if (!url1 || !url2 || !md5_1 || !md5_2 || qc_status == 'FILE_ISSUE') {
        println "Issue in metadata for ${run_accession}."
        return []
    }
   //${params.outDir}/$taxon_id/$run_accession/ 
    def pair1Path = "${params.outDir}/$taxon_id/$run_accession/${run_accession}_1.fastq.gz"
    def pair2Path = "${params.outDir}/$taxon_id/$run_accession/${run_accession}_2.fastq.gz"

//    def storedMd5 = getDataFromTable(run_accession,"data_file","md5").trim().split(';')
    def retryCount = 0
    def maxRetries = 3
    def md5Match = false

    while (!md5Match && retryCount < maxRetries) {
        // Download pair1
        
        """
        wget -qq -c -O ${pair1Path} ftp://${url1}
        """.execute().waitFor()
    
        """
        wget -qq -c -O ${pair2Path} ftp://${url2}
        """.execute().waitFor()
        // Create a shell script to download and calculate MD5
        /*
        def scriptContent = """
        wget -qq -c -O ${pair1Path} ftp://${url1}
        wget -qq -c -O ${pair2Path} ftp://${url2}
        """
        def scriptFile = new File("${params.outDir}/${taxon_id}/${run_accession}/download_script.sh")
        scriptFile.text = scriptContent
        scriptFile.setExecutable(true)
        */
        def md5Pair1 = DigestUtils.md5Hex(Files.newInputStream(file(pair1Path)))
        def md5Pair2 = DigestUtils.md5Hex(Files.newInputStream(file(pair2Path)))

        // Calculate MD5 checksums of downloaded files
        //def md5Pair1 = DigestUtils.md5Hex(file(pair1Path))
        //def md5Pair2 = DigestUtils.md5Hex(file(pair2Path))

        // Check if both MD5 checksums are present in stored MD5 checksums
        if ([md5_1, md5_2].containsAll([md5Pair1, md5Pair2])) {
            md5Match = true
            println "MD5 checksums match!"
        } else {
            println "MD5 checksums do not match! Retrying..."
            """
            rm ${pair1Path}
            """.execute().waitFor()
            """
            rm ${pair2Path}
            """.execute().waitFor()
            retryCount++
            Thread.sleep(1000) // Wait for 1 second before retrying
        }
    }

    if (!md5Match) {
        throw new RuntimeException("MD5 checksums do not match after $maxRetries retries!")
    }
    //fastqMetadata.add([taxon_id:taxon_id, gca:gca, run_accession:run_accession, file1:"${params.outDir}/${taxon_id}/${run_accession}/${run_accession}_1.fastq.gz", file2:"${params.outDir}/${taxon_id}/${run_accession}/${run_accession}_2.fastq.gz",query: "${params.outDir}/$taxon_id/$run_accession/${dataFileQuery}"])
    //"""
    //cp ${pair1Path} .
    //cp ${pair2Path} .
    //"""

    download_OUT=[]
    download_OUT.add([taxon_id:taxon_id, gca:gca, run_accession:run_accession, pair1:["${params.outDir}",taxon_id,run_accession,"${run_accession}_1.fastq.gz"].join("/"), pair2:["${params.outDir}",taxon_id,run_accession,"${run_accession}_2.fastq.gz"].join("/"),dataFileQuery:["${params.outDir}",taxon_id,run_accession,dataFileQuery].join("/")])
    """
    cp ${pair1Path} .
    cp ${pair2Path} .
    echo '${download_OUT.toString()}'
    """
    //return download_OUT
    /*
     // Perform the download and MD5 checksum verification
    """
    # Download the first file
    wget -O ${run_accession}_1.fastq.gz ${url1}
    if [ $? -ne 0 ]; then
        echo "Error downloading ${url1}"
        exit 1
    fi

    # Verify the MD5 checksum of the first file
    echo "${md5_1}  ${run_accession}_1.fastq.gz" | md5sum -c -
    if [ $? -ne 0 ]; then
        echo "MD5 checksum verification failed for ${run_accession}_1.fastq.gz"
        exit 1
    fi

    # Download the second file
    wget -O ${run_accession}_2.fastq.gz ${url2}
    if [ $? -ne 0 ]; then
        echo "Error downloading ${url2}"
        exit 1
    fi

    # Verify the MD5 checksum of the second file
    echo "${md5_2}  ${run_accession}_2.fastq.gz" | md5sum -c -
    if [ $? -ne 0 ]; then
        echo "MD5 checksum verification failed for ${run_accession}_2.fastq.gz"
        exit 1
    fi
    """
    */
}
