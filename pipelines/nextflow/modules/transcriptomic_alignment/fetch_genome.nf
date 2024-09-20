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

process FETCH_GENOME {
  tag "$taxon_id:$gca"
  label 'fetch_file'

  //storeDir "${params.outDir}/$taxon_id/$gca/ncbi_dataset/"
  afterScript "sleep $params.files_latency"  // Needed because of file system latency
  maxForks 10

  input:
  tuple val(taxon_id), val(gca), val(run_accession), val(pair1), val(pair2)
    

  output:
  tuple val(taxon_id), val(gca), val(run_accession), val(pair1), val(pair2), val("${params.outDir}/$taxon_id/$gca/ncbi_dataset/")
  
  script:
  """
  if [ ! -d "${params.outDir}/${taxon_id}/${gca}/ncbi_dataset/"  ]  ||  [ ! -n "\$(find "${params.outDir}/${taxon_id}/${gca}/ncbi_dataset/" -type f -name '*.fna' 2>/dev/null)" ]; then
    echo "Directory ncbi_dataset does not exist. Proceeding with download..."
    curl --retry 3  -X GET "${params.ncbiBaseUrl}/${gca}/download?include_annotation_type=GENOME_FASTA&hydrated=FULLY_HYDRATED" -H "Accept: application/zip" --output genome_file.zip
    unzip -j genome_file.zip
    mkdir -p ${params.outDir}/${taxon_id}/${gca}/ncbi_dataset && cp *.fna ${params.outDir}/${taxon_id}/${gca}/ncbi_dataset
  else
    echo "Directory ncbi_dataset already exists. Skipping download."
  fi
  """

}
