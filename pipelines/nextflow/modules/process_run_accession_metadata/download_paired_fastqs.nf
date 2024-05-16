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

import org.apache.commons.codec.digest.DigestUtils

process DOWNLOAD_PAIRED_FASTQS {
    label "default"
    tag "download ${run_accession} fastqs"
    maxForks 10
    publishDir "${params.outDir}/$taxon_id/$run_accession", mode: 'copy'

    input:
    tuple val(taxon_id), val(gca), val(run_accession)

    output:
    tuple (val(taxon_id), val(gca), val(run_accession), path(*_1.fastq.gz ), path(*_2.fastq.gz))

    when:
    len(fileUrls.trim().split(';')) == 2

    script:
    def fileUrls = getDataFileData(run_accession, "file_url")
    def (file1, file2) = fileUrls.trim().split(';')

    if (fileUrls.trim().split(';').size() != 2) {
        println "Expected two file URLs, but found ${fileUrls.size()}. Skipping the process."
        return
    }

    def pair1Path = "${publishDir}/${run_accession}_1.fastq.gz"
    def pair2Path = "${publishDir}/${run_accession}_2.fastq.gz"

    def storedMd5 = getDataFileData(run_accession, "md5").trim().split(';')
    def retryCount = 0
    def maxRetries = 3
    def md5Match = false

    while (!md5Match && retryCount < maxRetries) {
        // Download pair1
        "wget -qq -O ${pair1Path} ftp://${file1}".execute().waitFor()

        // Download pair2
        "wget -qq -O ${pair2Path} ftp://${file2}".execute().waitFor()

        // Calculate MD5 checksums of downloaded files
        def md5Pair1 = DigestUtils.md5Hex(new File(pair1Path))
        def md5Pair2 = DigestUtils.md5Hex(new File(pair2Path))

        // Check if both MD5 checksums are present in stored MD5 checksums
        if (storedMd5.containsAll([md5Pair1, md5Pair2])) {
            md5Match = true
            println "MD5 checksums match!"
        } else {
            println "MD5 checksums do not match! Retrying..."
            retryCount++
            Thread.sleep(1000) // Wait for 1 second before retrying
        }
    }

    if (!md5Match) {
        throw new RuntimeException("MD5 checksums do not match after $maxRetries retries!")
    }
}
