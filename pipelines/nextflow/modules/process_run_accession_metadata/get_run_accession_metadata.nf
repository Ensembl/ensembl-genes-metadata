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
// pyhton env?
process GET_RUN_ACCESSION_METADATA {

    label 'default'
    tag "$run_accession"

    input:
    tuple val(taxon_id), val(run_accession)

    output:
    file(joinPath(params.outDir, "${taxon_id}", "${run_accession}", "metadata.json")) into runAccessionsMetadataPath

    script:
    """
    // Construct the command based on whether last_date is provided
    def pythonScript = file("$projectDir/src/python/ensembl/genes/metadata/transcriptomic/get_metadata.py")
    def command = "python ${pythonScript} ${run_accession}"

    // Execute the Python script
    def process = ["python", pythonScript.toString(), run_accession].execute()
    process.waitFor()
    
    // Check if the script execution was successful
    if (process.exitValue() != 0) {
        throw new RuntimeException("Error executing Python script: ${pythonScript}")
    }

    // Get the output of the script
    def output = process.text.trim()

    // Verify if the returned path exists
    if (!new File(output).exists()) {
        throw new RuntimeException("The returned path does not exist: ${output}")
    }

    // Emit the path to the JSON file
    output
    """
}


