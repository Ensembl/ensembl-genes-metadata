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


// module description 
process STORE_METADATA {
    scratch false
    label 'default'
    tag 'store_metadata'
    
    input:
    path metadata2process
    val mysqlUpdate

    script:
    """
    log.info("Executing Python script to get metadata for run: $metadata2process")
    try {
    def cmd = "python write2db.py --file-path $run_accession --update $mysqlUpdate"
    def proc = cmd.execute()
    proc.waitFor()

    if (proc.exitValue() != 0) {
        throw new RuntimeException("Python script failed with exit code: ${proc.exitValue()}")
    }

    def output = proc.getText() // Capture the output of the process

    // Now you can use the 'output' variable as needed
    log.info("Output of the Python script: $output")

    // Splitting the output into individual parts
    def outputs = output.split('\n')

    // Checking if there's only one output or multiple outputs
    if (outputs.size() == 1) {
        // Only one output, pass it directly
        setMetaDataRecord(outputs[0])
    } else {
        // Multiple outputs, loop through and pass each one
        outputs.each { singleOutput ->
            setMetaDataRecord(singleOutput)
        }
    }
    } catch (Exception e) {
        log.error("Error executing Python script: ${e.message}")
        throw e // Rethrow the exception to halt the process
    }
    """
}