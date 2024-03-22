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


process FETCH_GENOME {

  label 'fetch_file'
  input:
  val gca
  val taxon_id
  val run_accession

  
  output:
  val taxon_id, emit:taxon_id
  val run_accession, emit:run_accession
  val(gca), emit:gca
  file(joinPath(params.outDir, "${taxon_id}", "${run_accession}", "ncbi_dataset", "data", "${gca}", "*.fna") ), emit: genome_file


  script:
  """
  curl -X GET "${params.ncbiBaseUrl}/${gca}/download?include_annotation_type=GENOME_FASTA&hydrated=FULLY_HYDRATED"  -H "Accept: application/zip" --output genome_file.zip
  unzip genome_file.zip
  
  //ncbi_dataset/data/GCA_963576655.1/GCA_963576655.1_icGasPoly1.1_genomic.fna 
  """
}