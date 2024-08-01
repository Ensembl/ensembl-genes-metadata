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
include { setMetaDataRecord } from './utils.nf'
process BUILD_QUERY {
    scratch false
    label 'python'
    tag "$run_accession"
    //maxForks 50

    input:
    tuple val(taxon_id), val(gca), val(run_accession), path(metadata2process)
    val mysqlUpdate
    
    output:
    tuple val(taxon_id), val(gca), val(run_accession), stdout
    //tuple val(taxon_id), val(gca), val(run_accession), val("${params.outDir}/${taxon_id}/${run_accession}/output_query.txt")
    //val(queryOutput)
//    stdout

    script:
    """
    chmod +x $projectDir/bin/write2db.py;
    $projectDir/bin/write2db.py --file_path ${metadata2process} --update ${mysqlUpdate} 
    """
    /*
     --output_dir ${params.outDir}/${taxon_id}/${run_accession}
     ${params.outDir}/${taxon_id}/${run_accession}/${metadata2process}
    output = Files.newInputStream(file("${params.outDir}/${taxon_id}/${run_accession}/output_query.txt"))
    queriesArray = output.text.trim().toString().replaceAll('; ', ' ').split(";")
    log.info("Queries Array: ${gca} ${metadata2process} ${queriesArray}")
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
    
    
    buildQueryPY=false
    queryOutput=[]
    queryOutput.add([taxon_id:taxon_id, gca:gca, run_accession:run_accession])
    log.info("Executing Python script to build query for run: $run_accession $metadata2process")
    outputFilePath = file("${params.outDir}/${taxon_id}/${run_accession}/output_query.txt")
    while (!outputFilePath.exists()){


    """
    chmod +x $projectDir/bin/write2db.py;
    """.execute().waitFor()
    """
    $projectDir/bin/write2db.py --file_path ${params.outDir}/${taxon_id}/${run_accession}/${metadata2process} --update ${mysqlUpdate} --output_dir ${params.outDir}/${taxon_id}/${run_accession} 
    """.execute().waitFor()
}
    output = Files.newInputStream(file("${params.outDir}/${taxon_id}/${run_accession}/output_query.txt")).text.trim()
    //queriesArray = output.tokenize(';')
    //println(queriesArray)
    // Split the output based on semicolons not within quotes or parentheses
    //queriesArray = output.split(/(?<![;\)\]])\s*;\s*(?![;\(\[])/).collect { it.trim() }
    //queriesArray = output.split(/(?<![;\)\]])\s*;\s*(?![;\(\[])/).collect { it.trim().replaceAll(/,\s+/, ', ') }
    queriesArray = output.toString().replaceAll('; ', ' ').split(";")
    //queriesArray = output.split(/(?<![\(\[]\);(?![^\(\[]*[\)\]])/).collect { it.trim() }
    //queriesArray =output.tokenize(/(?<!;)\s*;\s*(?!['"])/)
                    log.info("Queries Array: ${queriesArray}")
    // Determine if we have one or multiple queries
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
    log.info("ARRIVPO AL PRINT??")
    """
    echo '${queryOutput.toString()}'
    """
    */
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