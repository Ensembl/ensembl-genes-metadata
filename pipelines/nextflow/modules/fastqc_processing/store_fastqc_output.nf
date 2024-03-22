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
process STORE_FASTQC_OUTPUT {
    scratch false
    label 'default'
    tag 'store_fastqc_output'

    input:
    file fastqcMetadata

    script:
    """
    // Construct the command based on whether last_date is provided
    def pythonScript = file("$projectDir/src/python/ensembl/genes/metadata/database/JSON_to_db.py")
    def command = "python ${pythonScript} --file-path ${fastqcMetadata}"

    // Execute the Python script
    def process = command.execute()
    process.waitFor()
    
    // Check if the script execution was successful
    if (process.exitValue() != 0) {
        throw new RuntimeException("Error executing Python script: ${pythonScript}")
    }
    """
}