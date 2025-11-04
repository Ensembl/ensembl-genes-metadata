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

process SET_DATE {
    output:
    stdout

    script:
    if(params.date && !params.full_screen)
        """
            echo ${params.date}
        """
    else if(params.full_screen && !params.date)
        """
            ${params.metadata_db.host} ${params.metadata_db.database}  -NB -e "SELECT DATE_FORMAT(date_value, '%m/%d/%Y') from update_date WHERE update_type = 'full_screen';"
        """
    else if(!params.date && !params.full_screen)
        """
            ${params.metadata_db.host} ${params.metadata_db.database}  -NB -e "SELECT DATE_FORMAT(DATE_SUB(date_value, INTERVAL 1 DAY), '%m/%d/%Y') from update_date WHERE update_type = 'regular_update';"
        """ 
    else
        error "Invalid parameters to set up date"

}