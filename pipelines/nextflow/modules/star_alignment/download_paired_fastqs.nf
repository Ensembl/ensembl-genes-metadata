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


process DOWNLOAD_PAIRED_FASTQS {
    label "default"
    tag "${taxon_id}:${run_accession}"
    maxForks 25
    //storeDir "${params.outDir}/$taxon_id/$run_accession"
    afterScript "sleep $params.files_latency"  // Needed because of file system latency
    
    input:
    tuple val(taxon_id), val(gca), val(run_accession), val(url1), val(url2), val(md5_1), val(md5_2), val(genomeDir)

    output:
    tuple val(taxon_id), val(gca), val(run_accession), val("${params.outDir}/$taxon_id/$run_accession/${run_accession}_1.fastq.gz"), val("${params.outDir}/$taxon_id/$run_accession/${run_accession}_2.fastq.gz" ), val(genomeDir)
    

    script: 
    if(!url1 || !url2 || !md5_1 || !md5_2){
            log.error("Metadata corrupted")
    }else{
        def pair1Path = "${params.outDir}/$taxon_id/$run_accession/${run_accession}_1.fastq.gz"
        def pair2Path = "${params.outDir}/$taxon_id/$run_accession/${run_accession}_2.fastq.gz"
        def retryCount = 0
        def maxRetries = 3
        def md5Match = false
        def filesValid = false
        def dirPath = file("${params.outDir}/${taxon_id}/${run_accession}")
        // List and filter .gz files
        def gzFiles = dirPath.listFiles().findAll { it.name.endsWith('.gz') }

        log.info("gzFiles ${gzFiles}")
        if (gzFiles.size() !=2 ){
            filesValid = true
        
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
            log.info("MD5 checksums do not match after $maxRetries retries!")
        } 
        }
    }
    
    """
    echo "Download completed"
    """
}
