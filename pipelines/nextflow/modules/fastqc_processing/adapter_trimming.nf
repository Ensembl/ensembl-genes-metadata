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
process ADAPTER_TRIMMING {
    scratch true
    label 'default'
    tag 'trimming '

    input:
    tuple val(taxon_id), val(run_accession), val(pairedFastqFiles)
    set pair1, pair2 from pairedFastqFiles

    output:
    tuple (taxon_id, run_accession, trimmedFastqFiles) into runTrimmedFastqs
    //no need to emit the path because the subsampled files will be in the run_accession dir _1_1

    script:
    """
    // Construct the command based on whether last_date is provided
    def pythonScript = file("$projectDir/src/python/ensembl/genes/ ")
    def command = "python ${pythonScript} "

    // Execute the Python script
    def process = command.execute()
    process.waitFor()

    // Check if the script execution was successful
    if (process.exitValue() != 0) {
        throw new RuntimeException("Error executing Python script: ${pythonScript}")
    }

    // Read the output of the Python script
    def output = []
    process.inputStream.eachLine { line ->
            output.add(line.trim())
        }

    // Define the sampledFastqFiles array
    def trimmedFastqFiles = []

    // Add the captured paths from the output array to sampledFastqFiles
    output.each { path ->
        trimmedFastqFiles.add(new File(path))
    }

    // Check if the correct number of files were captured
    if (trimmedFastqFiles.size() != 2) {
        throw new RuntimeException("Expected two paths from Python script, but received ${sampledFastqFiles.size()}")
    }

    // Emit the tuple
    emit(taxon_id, run_accession, trimmedFastqFiles)
    """
}