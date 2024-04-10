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

process GET_RUN_ACCESSION {

    label 'default'
    tag "$taxon_id"

    input:
    val taxon_id
    val dateQuery = null // Optional input parameter for date query (format: 'YYYY-MM-DD')

    output:
    file(joinPath(params.outDir, "${taxon_id}", "run_accessions.txt")) into runAccessionsPath

    script:
    """
    def taxonQuery = "tax_eq(${taxon_id})"
    def instrumentQuery = "instrument_platform=ILLUMINA"
    def layoutQuery = "library_layout=PAIRED"
    def sourceQuery = "library_source=TRANSCRIPTOMIC"
    def defaultDateQuery = "first_created >= '2019-01-01'" // Default date query if dateQuery is not provided

    // Use provided dateQuery if available, otherwise use default
    def usedDateQuery = dateQuery ? "first_created >= '${dateQuery}'" : defaultDateQuery

    def query = [taxonQuery, instrumentQuery, layoutQuery, sourceQuery, usedDateQuery].join(" AND ")
    def encodedQuery = URLEncoder.encode(query, 'UTF-8')

    def searchUrl = "https://www.ebi.ac.uk/ena/portal/api/search?display=report&query=${encodedQuery}&domain=read&result=read_run&fields=run_accession"

    def response = new URL(searchUrl).text
    def runAccessions = response.split("\\n").drop(1).join("\\n") // Remove header

    new File('run_accessions.csv').write(runAccessions)
    """
}
