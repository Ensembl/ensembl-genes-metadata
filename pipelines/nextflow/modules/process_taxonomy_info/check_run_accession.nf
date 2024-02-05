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
process CHECK_RUN_ACCESSION {
    scratch false
    label 'default'
    tag "$run_accession"
    
    input:
    val taxon_id
    
    val transcriptomic_dbname, 
    val  transcriptomic_host,
    val  transcriptomic_port,   
    val  transcriptomic_user,
    val  transcriptomic_password,
    val run_accession_list

    output:
    val filtered_run_accessions, emit : filtered_run_accessions

    script:
    """
    """
    emit 
}