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

include { getRunTable } from '../modules/utils.nf'
include { getDataFromTable } from '../modules/utils.nf'

process EXTRACT_UNIQUELY_MAPPED_READS_PERCENTAGE {
    label 'python'
    tag "percentage mapped in  $run_accession"
    publishDir "${params.outDir}/$taxon_id/$run_accession/star/", mode: 'copy'
    afterScript "sleep $params.files_latency"  // Needed because of file system latency

    input:
    tuple val(taxon_id), val(gca), val(run_accession), val(starOutputFile)


    output:
    tuple val(taxon_id), val(gca), val(run_accession)
    path("insert_into_align.json")

    script:
    def run_id = getDataFromTable("run_id", "run", "run_accession", run_accession)[0].run_id
    //def run_Id = getRunTable(run_accession, 'run_id')
    //def star_dir = "${params.outDir}/${taxon_id}/${run_accession}/star/"
    def star_dir = new File(starOutputFile).parent
    """
    chmod +x parse_star_output.py
    parse_star_output.py --file_path ${starOutputFile} --output_dir ${star_dir}\
    --extra_parameters "{'run_accession': '${run_accession}', 'assembly_accession': '${gca}', 'run_id': '${run_id}'}" \
    --uniquely_mapped_reads_percentage --percentage_reads_mapped_to_multiple_loci --percentage_reads_unmapped_too_short"
    """
}
