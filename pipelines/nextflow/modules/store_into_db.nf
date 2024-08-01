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
import java.nio.file.Files
import java.nio.file.Paths
import groovy.io.FileType
 import java.util.concurrent.TimeUnit
include { setMetaDataRecord } from './utils.nf'
process STORE_INTO_DB {
    scratch false
    label 'default'
    tag "$run_accession"
    //maxForks 20

    input:
    tuple val(taxon_id), val(gca), val(run_accession), val(query)
    
    output:
    tuple val(taxon_id), val(gca), val(run_accession)

    script:
//    int maxRetries = 5
  //  int retryCount = 0
    //while (retryCount < maxRetries) {
      //  if (file(query).exists()) {
        //    break
       // }
       // println "File ${query} does not exist yet. Retrying in 5 seconds..."
       // TimeUnit.SECONDS.sleep(5)
        //retryCount++
    //}
    //if (retryCount == maxRetries) {
      //  println "File ${query} was not found after ${maxRetries} attempts. Exiting."
       // exit 1
    //}
    //output = new File(query).text.trim()
    queriesArray = query.replaceAll('; ', ' ').split(";")
    log.info("Queries Array:  ${run_accession}  ${queriesArray} ")
    if (queriesArray.size() == 1) {
    setMetaDataRecord(queriesArray[0].trim().toString())
    } else if (queriesArray.size() > 1) {
    queriesArray.each { query ->
            // Trim the query to remove any leading/trailing whitespace
            query = query.trim()
            // Check if the query is not empty
            if (query) {
                //log.info("queriesArray[${index}] ${query};")
                setMetaDataRecord(query.toString())
            }
        }
    }

    """
    echo "${queriesArray.toString().trim()}"
    """
    /*
    // Loop through the queriesArray and log information
        queriesArray.each { query ->
                log.info("queriesArray ${query};")
                // Trim the query to remove leading/trailing whitespace (if necessary)
                    query = query
                                if (query) {

                        setMetaDataRecord(query + ";")
                            }
                            }
                            

    def outputDirPath = "${params.outDir}/${taxon_id}/${run_accession}"
    def outputDir = new File(outputDirPath)
    def pattern = ~/.*output_query\.txt/

    if (outputDir.exists() && outputDir.isDirectory()) {
    def matchingFiles = []
    outputDir.eachFileMatch(FileType.FILES, pattern) { file ->
        matchingFiles << file
    }

    if (!matchingFiles.isEmpty()) {
        matchingFiles.each { file ->
            log.info("Found file: ${file.absolutePath}")

            // Use Files.newInputStream to read the file
            def path = Paths.get(file.absolutePath)
            def output = ""
            Files.newInputStream(path).withCloseable { inputStream ->
                output = new BufferedReader(new InputStreamReader(inputStream)).lines().collect(Collectors.joining("\n")).trim()
            }
            
            def queriesArray = output.tokenize(';')
            println(queriesArray)
            // Loop through the queriesArray and log information
            queriesArray.each { query ->
                log.info("queriesArray ${query};")
                // Trim the query to remove leading/trailing whitespace (if necessary)
                query = query.trim()
                setMetaDataRecord(query.toString() + ";")
             }
        }
        }
    
    output = Files.newInputStream(file("${params.outDir}/${taxon_id}/${run_accession}/output_query.txt")).text.trim()
    queriesArray = output.tokenize(';')
    println(queriesArray)
    // Loop through the queriesArray and log information
        queriesArray.each { query ->
                log.info("queriesArray ${query};")
                // Trim the query to remove leading/trailing whitespace (if necessary)
                    query = query.trim()
                        setMetaDataRecord(query.toString() + ";")
                            }
                            

    for (int index = 0; index < queriesArray.size(); index++) {
        query = queriesArray[index]
            log.info("queriesArray[${index}] ${query};")
                setMetaDataRecord(query.toString() + ";")
                }
    
    
    #def queriesArray = output.toString().split(";")
    queriesArray.eachWithIndex { query, index ->
         // Trim the query to remove any leading/trailing whitespace
         query = query.trim()
         // Check if the query is not empty
         if (query) {
             log.info("queriesArray ${query};")
             setMetaDataRecord(query.toString() + ";")                                                                     

                   }

                }
    */            
}
