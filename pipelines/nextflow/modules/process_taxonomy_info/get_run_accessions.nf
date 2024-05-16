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
// Define a class or data structure to hold the accumulated lines
class ResponseData {
    List<String> lines = []
}

include { setMetaRecord } from '../utils.nf'
process GET_RUN_ACCESSIONS {
    label 'default'
    tag "$taxon_id:$gca"
    storeDir "${params.outDir}/$taxon_id/"
    publishDir "${params.outDir}/$taxon_id", mode: 'copy'
    maxForks 10
    afterScript "sleep $params.files_latency"  // Needed because of file system latency

    input:
    tuple val(taxon_id), val(gca), val(lastCheckedDate)

    output:
    //val(runAccessionChannel), emit: runAccessionData
    val(runAccessionList)
    path("run_accession_list.txt")
    
    script:
  //  println("GET RUN ACCESSION")
    def taxonQuery = "tax_eq(${taxon_id})"
    def instrumentQuery = "instrument_platform=ILLUMINA"
    def layoutQuery = "library_layout=PAIRED"
    def sourceQuery = "library_source=TRANSCRIPTOMIC"
    // Use provided dateQuery if available, otherwise use default
    def usedDateQuery = "first_created >='${lastCheckedDate.trim()}'"

    def searchUrl = "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=${taxonQuery}%20AND%20instrument_platform=ILLUMINA%20AND%20library_layout=PAIRED%20AND%20library_source=TRANSCRIPTOMIC%20AND%20first_created%3E=${lastCheckedDate.trim()}&domain=read&fields=run_accession"

    //try {
    // Open a connection to the URL
    URL url = new URL(searchUrl)
    URLConnection connection = url.openConnection()

    // Set the connection and read timeouts (in milliseconds)
    connection.setConnectTimeout(5000) // 5 seconds
    connection.setReadTimeout(10000)   // 10 seconds

    // Get the input stream
    BufferedReader reader = new BufferedReader(new InputStreamReader(connection.getInputStream()))
    
    // Create a list to hold the run accession data
    //List<Map<String, String>> runAccessionList = []
    runAccessionList = []
    // Create a string to accumulate the response content
    StringBuilder responseContent = new StringBuilder()

    // Process the response line by line
    String line
    while ((line = reader.readLine()) != null){
        // Append each line to the response content

        responseContent.append(line).append('\n')
        
        if (!line.startsWith("run_accession")) { // Skip header line
            runAccessionList.add([taxon_id: taxon_id, gca: gca, run_accession: line])
    //    println("Response7777: $line")
        }
      //  println("Response: $runAccessionList")
    }
    reader.close()

//    responseContent.toString().split('\n').each { row ->
  //  if (!row.startsWith("run_accession")) {
        //if (row !="run_accession") { // Skip header line
    //        runAccessionList.add(['taxon_id': taxon_id, 'gca': gca, 'run_accession': row])
      //      println("Response: $runAccessionList")
       // }
       // }
    // Save the response content to a text file
    //def responseFile = new File('response.txt')
    //responseFile.text = responseContent.toString()
//    println(runAccessionList)
    // Create a channel from the list of run accessions
//    runAccessionChannel = Channel.from(runAccessionList)
    // Create a copy of runAccessionList
  //  def runAccessionListCopy = runAccessionList.toList()

    // Create a channel from the copied list
//    runAccessionChannel = Channel.from(runAccessionListCopy)
//    runAccessionChannel.view { item -> println(item) }
    runAccessionToFile='run_accession_list.txt'
    """
    echo '${responseContent.toString()}' > $runAccessionToFile
    """
//echo '${responseContent.toString()}' > $runAccessionToFile
//    def runAccessionList = []
    //def file = new File('run_accessions.csv')
    //file.withWriter { writer ->
//        responseData.lines.each { lineData ->
        //file.eachLine { line -> writer.writeLine([taxon_id: ${taxon_id}, gca: ${gca}, run_accession: lineData])}
//       runAccessionList.add([taxon_id: taxon_id, gca: gca, run_accession: lineData])
     //   writer.writeLine(lineData)
  //      println("Response: $lineData")
    //}
   // }

    //runAccessionData = Channel.from(responseData.lines)
      //                 .map { lineData -> [taxon_id: taxon_id, gca: gca, run_accession: lineData]}

   //setMetaRecord(taxon_id,'update')
//runAccessionChannel.flatten().view { d -> "AAAA Taxon ID: ${d.taxon_id}, GCA: ${d.gca}, last date: ${d.run_accession}"} 
//    } catch (Exception e) {
    // Handle any exceptions
  //  println("Error: ${e.message}")
   // }
}
