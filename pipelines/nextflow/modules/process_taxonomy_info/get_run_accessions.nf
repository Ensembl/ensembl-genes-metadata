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
@Grab('org.codehaus.groovy:groovy-all:2.2.2')
import java.net.URLEncoder
import java.net.URLConnection
import java.net.URL
import java.net.URLConnection
import java.io.BufferedReader
import java.io.InputStreamReader
import java.io.File
import groovy.json.JsonOutput

process GET_RUN_ACCESSIONS {
    label 'default'
    tag "$taxon_id:$gca"
    storeDir "${params.outDir}/$taxon_id/"
    publishDir "${params.outDir}/$taxon_id", mode: 'copy'
    maxForks 10
    afterScript "sleep $params.files_latency"  // Needed because of file system latency

    input:
    tuple val(taxon_id), val(gca), val(run_accession_batch), val(lastCheckedDate)

    output:
    val(runAccessionList)
    path("run_accession_list.txt")
    
    script:
    runAccessionList = []    
    runAccessionToFile='run_accession_list.txt' 
    //def fileBatch=new File(run_accession_batch)
    //if (fileBatch.exists()){
    if (run_accession_batch && file(run_accession_batch).exists()){
        fileBatch.eachLine { line ->
                runAccessionList.add([taxon_id: taxon_id, gca: gca, run_accession: line])
    }   
    """
    cat ${run_accession_batch.toString()} > $runAccessionToFile
    """
    } else {
    def taxonQuery = "tax_eq(${taxon_id})"
    def instrumentQuery = "instrument_platform=ILLUMINA"
    def layoutQuery = "library_layout=PAIRED"
    def sourceQuery = "library_source=TRANSCRIPTOMIC"
    def usedDateQuery = "first_created >='${lastCheckedDate.trim()}'"

    def searchUrl = "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=${taxonQuery}%20AND%20instrument_platform=ILLUMINA%20AND%20library_layout=PAIRED%20AND%20library_source=TRANSCRIPTOMIC%20AND%20first_created%3E=${lastCheckedDate.trim()}&domain=read&fields=run_accession"

    // Open a connection to the URL
    URL url = new URL(searchUrl)
    URLConnection connection = url.openConnection()

    // Set the connection and read timeouts (in milliseconds)
    connection.setConnectTimeout(10000) // 5 seconds
    connection.setReadTimeout(10000)   // 10 seconds

    // Get the input stream
    BufferedReader reader = new BufferedReader(new InputStreamReader(connection.getInputStream()))
    
    // Create a list to hold the run accession data
    //runAccessionList = []
    // Create a string to accumulate the response content
    StringBuilder responseContent = new StringBuilder()

    // Process the response line by line
    String line
    while ((line = reader.readLine()) != null){
        // Append each line to the response content
        responseContent.append(line).append('\n')
        
        if (!line.startsWith("run_accession")) { // Skip header line
            runAccessionList.add([taxon_id: taxon_id, gca: gca, run_accession: line])
        }
    }
    reader.close()

    log.info("BATCH FILE ${run_accession_batch}")
    """ 
    echo '${responseContent.toString()}' > $runAccessionToFile
    """
}
}
