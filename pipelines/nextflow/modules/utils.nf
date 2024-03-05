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
import groovy.sql.Sql
import java.time.LocalDateTime
import java.time.format.DateTimeFormatter

def checkTaxonomy(String jdbcUrl, String username, String password, String taxonId) {
    def sql = Sql.newInstance(jdbcUrl, username, password)
    
    try {
        def query = "SELECT * FROM meta  WHERE taxon_id = '${taxonId}'"
        def result = sql.rows(query)
        return result.size() > 0
    } finally {
        sql.close()
    }
}

def getLastCheckDate(String jdbcUrl, String username, String password, String taxonId) {
    def sql = Sql.newInstance(jdbcUrl, username, password)
    def lastCheckDate = null

    try {
        def query = "SELECT last_check FROM meta WHERE taxon_id = '${taxonId}'"
        def result = sql.rows(query)

        if (result.size() > 0) {
            // Assuming 'last_check' is a date-like column
            // Adjust the date format pattern based on the actual format in your database
            def dateFormat = new SimpleDateFormat("yyyy-MM-dd") // Adjust the format if needed
            lastCheckDate = dateFormat.parse(result[0].last_check)
        }
    } finally {
        sql.close()
    }

    return lastCheckDate
}

def insertMetaRecord(String jdbcUrl, String username, String password, String taxonId) {
    def sql = Sql.newInstance(jdbcUrl, username, password)

    try {
        // Get the current date and time
        def currentDate = LocalDateTime.now()
        def dateFormatter = DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss")
        def formattedDate = currentDate.format(dateFormatter)

        // Execute the SQL INSERT statement
        def insertQuery = "INSERT INTO meta (taxon_id, last_check) VALUES ('${taxonId}', '${formattedDate}')"
        sql.executeUpdate(insertQuery, 'meta_id')

    } finally {
        sql.close()
    }

}

def build_ncbi_path(gca, assembly_name) {
    final gca_splitted = gca.replaceAll("_","").tokenize(".")[0].split("(?<=\\G.{3})").join('/')
    return  'https://ftp.ncbi.nlm.nih.gov/genomes/all'  + '/' + gca_splitted + '/' + "$gca" +'_' + assembly_name.replaceAll(" ","_") + '/' + "$gca" + '_' + assembly_name.replaceAll(" ","_") + '_genomic.fna.gz'
}


def concatString(string1, string2, string3){
    return string1 + '_'+string2 + '_'+string3
}

def getDataset(busco_lineage){
    String fileContents = new File(busco_lineage).text
    return fileContents
}