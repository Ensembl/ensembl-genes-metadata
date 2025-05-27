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

nextflow.enable.dsl=2

// Load Plugins
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Subworkflows
include {ALIGNMENT_PIPELINE} from "./pipelines/nextflow/workflows/alignment_pipeline.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow  {

    main:
    log.info "Pipeline started at: ${new Date().format('dd-MM-yyyy HH:mm:ss')}"
    

    // Validate input parameters
    validateParameters()
    // Print summary of supplied parameters
    log.info paramsSummaryLog(workflow)
    // Print summary of supplied parameters
    ALIGNMENT_PIPELINE(params.csvFile, params.bam2cram, params.mergeTissue, params.stranded, params.bam2bigWig)
    

workflow.onComplete {
    log.info "Pipeline completed at: ${new Date().format('dd-MM-yyyy HH:mm:ss')}"
    log.info  "Time to complete workflow execution: ${workflow.duration}"
    log.info  "Execution status: ${workflow.success ? 'Successful' : 'Failed'}"
    if (params.cleanOutputDir) {    
    try {
        def outDir = java.nio.file.Paths.get(params.outDir)

        java.nio.file.Files.newDirectoryStream(outDir, "*").each { path ->
        if (!path.toString().endsWith(".gz")) {
            deleteRecursively(path)
            }
        }
        log.info "Cleaning process completed successfully."
    } catch (Exception e) {
        log.error "Exception occurred while executing cleaning command: ${e.message}", e
    }
    }
}
workflow.onError {
    println "Error... Pipeline execution stopped with the following message: $workflow.errorMessage"
}
}
def deleteRecursively(Path path) {
    if (java.nio.file.Files.isDirectory(path)) {
        java.nio.file.Files.newDirectoryStream(path).each { subPath ->
            deleteRecursively(subPath)
        }
    }
    java.nio.file.Files.delete(path)
}

