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
// pyhton env?
process GET_RUN_ACCESSION_METADATA {

    label 'default'
    tag "$run_accession"
    publishDir "${params.outDir}/$taxon_id/$run_accession", mode: 'copy'
    maxForks 10
    input:
    tuple val(taxon_id), val(gca), val(run_accession)

    output:
    tuple val(taxon_id), val(gca), val(run_accession)
    path("insert_into_run.json")
    path("insert_into_study.json")
    path("insert_into_data_file.json")

    script:
    script:
    """
    run_accession = '${run_accession}' // Access run_accession from input channel or parameters

    log.info("Executing Python script to get metadata for run: $run_accession")
    try {
        def cmd = "python get_metadata.py --run $run_accession"
        def proc = cmd.execute()
        proc.waitFor()

        if (proc.exitValue() != 0) {
            throw new RuntimeException("Python script failed with exit code: ${proc.exitValue()}")
        }
    } catch (Exception e) {
        log.error("Error executing Python script: ${e.message}")
        throw e // Rethrow the exception to halt the process
    }
    """
}


