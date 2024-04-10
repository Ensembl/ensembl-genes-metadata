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


process EXTRACT_UNIQUELY_MAPPED_READS_PERCENTAGE {
    
    input:
    file starOutputFile 
    val run_accession
    val gca

    output:
    file "${params.outDir}/${taxon_id}/${run_accession}/star_alignment/metadata.json", emit: star_metadata


    script:
    """
    def run_id=getRunId(params.jdbcUrl, params.transcriptomic_dbuser, params.transcriptomic_dbpassword, $run_accession)
    def star_dir = new File(input_file).parent

    // Construct the command based on whether last_date is provided
    def pythonScript = file("$projectDir/ src/python/ensembl/genes/metadata/transcriptomic/parse_star_output.py")
    def command = "python ${pythonScript} --file_path $starOutputFile --output_dir $star_dir \
    --extra_parameters "{'run_accession': '$run_accession', 'assembly_accession': '$gca', 'run_id': '$run_id'}" \
    --uniquely_mapped_reads_percentage --percentage_reads_mapped_to_multiple_loci --percentage_reads_unmapped_too_short"

    // Execute the Python script
    def process = command.execute()
    process.waitFor()
    
    // Check if the script execution was successful
    if (process.exitValue() != 0) {
        throw new RuntimeException("Error executing Python script: ${pythonScript}")
    }
    """
}
