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

import groovy.json.JsonSlurper
import org.apache.commons.codec.digest.DigestUtils
import java.nio.file.Files
import java.nio.file.Path

include { updateTable } from '../utils.nf'
include { getDataFromTable } from '../utils.nf'

process DOWNLOAD_PAIRED_FASTQS {
    label "default"
    tag "${taxon_id}:${run_accession}"
    maxForks 25
    //storeDir "${params.outDir}/$taxon_id/$run_accession"
    afterScript "sleep $params.files_latency"  // Needed because of file system latency
    
    input:
    tuple val(taxon_id), val(gca), val(run_accession)
    //path(dataFileQuery)

    output:
    //val(download_OUT)
    tuple val(taxon_id), val(gca), val(run_accession), val("${params.outDir}/$taxon_id/$run_accession/${run_accession}_1.fastq.gz"), val("${params.outDir}/$taxon_id/$run_accession/${run_accession}_2.fastq.gz" ), val("${params.outDir}/$taxon_id/$run_accession/insert_into_data_file.json")
    
    //!file("${params.outDir}/${taxon_id}/${run_accession}/${run_accession}_1.fastq.gz}").exists() ||
    //!file("${params.outDir}/${taxon_id}/${run_accession}/${run_accession}_2.fastq.gz").exists() 
   script: 
    // Parse the JSON data
    def parsedJson = new JsonSlurper().parse(file("${params.outDir}/$taxon_id/$run_accession/insert_into_data_file.json"))
    def dataFiles = parsedJson.data_files
    //def download_OUT = [[taxon_id:taxon_id, gca:gca, run_accession:run_accession, pair1:["${params.outDir}",taxon_id,run_accession,"${run_accession}_1.fastq.gz"].join("/"), pair2:["${params.outDir}",taxon_id,run_accession,"${run_accession}_2.fastq.gz"].join("/"),dataFileQuery:["${params.outDir}",taxon_id,run_accession,dataFileQuery].join("/")]]

    // Check if there are exactly two data files
    //if (dataFiles.size() != 2) {
    //    return """
    //    echo '${download_OUT.toString()}'
    //    """
    //}
    if (dataFiles.size() == 2){
    def file1 = dataFiles[0]
    def file2 = dataFiles[1]
    println(file1)
    println(file2)
    // Extract URLs and MD5 checksums
    def url1 = file1.file_url
    def md5_1 = file1.md5
    def url2 = file2.file_url
    def md5_2 = file2.md5
    if(!url1 || !url2 || !md5_1 || !md5_2){
            log.error("Metadata corrupted")
    }else{
    def pair1Path = "${params.outDir}/$taxon_id/$run_accession/${run_accession}_1.fastq.gz"
    def pair2Path = "${params.outDir}/$taxon_id/$run_accession/${run_accession}_2.fastq.gz"
    def pair1SubPath = "${params.outDir}/$taxon_id/$run_accession/${run_accession}_1.fastq.gz.sub"
    def pair2SubPath = "${params.outDir}/$taxon_id/$run_accession/${run_accession}_2.fastq.gz.sub"
    def retryCount = 0
    def maxRetries = 3
    def md5Match = false
    def filesValid = false
    //def qc_status = getDataFromTable("qc_status", "run", "run_accession", run_accession)[0].qc_status
    def dirPath = file("${params.outDir}/${taxon_id}/${run_accession}")
    // List and filter .gz files
    def gzFiles = dirPath.listFiles().findAll { it.name.endsWith('.gz') }
    // List and filter .gz.sub files
    def subFiles = dirPath.listFiles().findAll { it.name.endsWith('.gz.sub') }
    log.info("subFiles ${subFiles} gzFiles ${gzFiles}")
    if (gzFiles.size() !=2 || subFiles.size() !=2){
        filesValid = true
      //  return """
      //  echo '${download_OUT.toString()}'
      //  """
    //}
    // Check for file issues and QC status
    //if (!url1 || !url2 || !md5_1 || !md5_2) { // || qc_status == 'FILE_ISSUE') {
    //    return 
    //    """
    //    echo '${download_OUT.toString()}'
    //    """
    //}

    while (!md5Match && retryCount < maxRetries) {
        """
        wget -qq -c -O ${pair1Path} ftp://${url1}
        """.execute().waitFor()
    
        """
        wget -qq -c -O ${pair2Path} ftp://${url2}
        """.execute().waitFor()

        // Calculate MD5 checksums of downloaded files
        md5Pair1 = DigestUtils.md5Hex(Files.newInputStream(file(pair1Path)))
        md5Pair2 = DigestUtils.md5Hex(Files.newInputStream(file(pair2Path)))

        // Check if both MD5 checksums are present in stored MD5 checksums
        if ([md5_1, md5_2].containsAll([md5Pair1, md5Pair2])) {
            md5Match = true
            log.info("MD5 checksums match!")
        } else {
            log.info("MD5 checksums do not match! Retrying...")
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
        //throw new RuntimeException("MD5 checksums do not match after $maxRetries retries!")
        log.info("MD5 checksums do not match after $maxRetries retries!")
        updateTable("run_accession", run_accession, "run", "qc_status", "DOWNLOAD_FAILED")
        //return """
        //echo '${download_OUT.toString()}'
        //"""
    } else {
    //download_OUT=[]
    //download_OUT.add([taxon_id:taxon_id, gca:gca, run_accession:run_accession, pair1:["${params.outDir}",taxon_id,run_accession,"${run_accession}_1.fastq.gz"].join("/"), pair2:["${params.outDir}",taxon_id,run_accession,"${run_accession}_2.fastq.gz"].join("/"),dataFileQuery:["${params.outDir}",taxon_id,run_accession,dataFileQuery].join("/")])
    """
    cp ${pair1Path} .
    cp ${pair2Path} .
    """
       }
    }
    }
    }
    """
    echo "Download completed"
    """
}
