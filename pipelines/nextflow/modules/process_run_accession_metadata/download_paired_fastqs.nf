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

process DOWNLOAD_PAIRED_FASTQ {
    label "default"
    tag "download ${run_accession} fastqs"
    // Define input parameters
    input:
    val taxon_id
    val run_accession

    // Define output channel
    output:
    set pair1, pair2 into pairedFastqFiles
    

    script:
    """
    def srr = run_accession.take(6)
    def first = srr.take(6)
    def second = '00' + srr[-1]
    def third = '0' + srr[-2..-1].join()

    # Extract the SRR prefix from the run accession
    srr_prefix=\$(echo ${run_accession} | cut -c1-6)
    # Extract the last digit from the run accession
    last_digit=\$(echo ${run_accession} | rev | cut -c1 | rev)

    # Define the base URL for the FTP server
    
    def pair1 = "wget -qq ${params.ftpBaseUrl}/${first}/${second}/${srr}/${run_accession}_1.fastq.gz -P ${path}".execute().waitFor()
    def pair2 = "wget -qq ${params.ftpBaseUrl}/${first}/${second}/${srr}/${run_accession}_2.fastq.gz -P ${path}".execute().waitFor()

    if (pair1 != 0) {
        // Retry with different URL pattern
        pair1 = "wget -qq ${params.ftpBaseUrl}/${first}/${srr}/${run_accession}_1.fastq.gz -P ${path}".execute().waitFor()
        pair2 = "wget -qq ${params.ftpBaseUrl}/${first}/${srr}/${run_accession}_2.fastq.gz -P ${path}".execute().waitFor()

        if (pair1 != 0) {
            // Retry with another URL pattern
            pair1 = "wget -qq ${params.ftpBaseUrl}/${first}/${third}/${srr}/${run_accession}_1.fastq.gz -P ${path}".execute().waitFor()
            pair2 = "wget -qq ${params.ftpBaseUrl}/${first}/${third}/${srr}/${run_accession}_2.fastq.gz -P ${path}".execute().waitFor()

            if (res != 0) {
                // If all attempts fail, throw an error
                throw new RuntimeException("Could not download ${run_accession}_1.fastq.gz")
            }
        }
    }

    # Download the paired FASTQ files
    wget -qq -P joinPath(params.outDir, "${taxon_id}", "${run_accession}") \${pair1}
    wget -qq -P joinPath(params.outDir, "${taxon_id}", "${run_accession}") \${pair2}
  
    // Construct file paths for pair1 and pair2
    def pair1Path = joinPath(params.outDir, "${taxon_id}", "${run_accession}", "${run_accession}_1.fastq.gz")
    def pair2Path = joinPath(params.outDir, "${taxon_id}", "${run_accession}", "${run_accession}_2.fastq.gz")

    // Emit both paths
    emit pair1: file(pair1Path), pair2: file(pair2Path)
    """
}

