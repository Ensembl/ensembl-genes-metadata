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
import groovy.sql.Sql
import java.time.LocalDateTime
import java.time.format.DateTimeFormatter

def checkTaxonomy(String taxonId) {
    def sql 
    def driver = 'com.mysql.cj.jdbc.Driver'
    def jdbcUrl = "jdbc:mysql://${params.transcriptomic_dbhost}:${params.transcriptomic_dbport}/${params.transcriptomic_dbname}"
    sql = Sql.newInstance(jdbcUrl, params.transcriptomic_dbuser,params.transcriptomic_dbpassword,driver)
    
    try {
        def query = "SELECT * FROM meta  WHERE taxon_id = ? "
        def result = sql.rows(query,[taxonId])
        return result.size() > 0
    } catch (Exception ex) {
        ex.printStackTrace()}
    finally {
        sql.close()
    }
}

def getLastCheckedDate(String taxonId) {
    def sql 
    def driver = 'com.mysql.cj.jdbc.Driver'
    def jdbcUrl = "jdbc:mysql://${params.transcriptomic_dbhost}:${params.transcriptomic_dbport}/${params.transcriptomic_dbname}"
    sql = Sql.newInstance(jdbcUrl, params.transcriptomic_dbuser,params.transcriptomic_dbpassword,driver)
    def lastCheckedDate = null

    try {
        def query = "SELECT last_check FROM meta WHERE taxon_id = ? "
        def result = sql.rows(query, [taxonId])
        if (result.size() > 0) {
            // Assuming 'last_check' is a date-like column
            // Adjust the date format pattern based on the actual format in your database
            def dateFormat = new SimpleDateFormat("yyyy-MM-dd") // Adjust the format if needed
            lastCheckedDate = dateFormat.parse(result[0].last_check)
        }
    } catch (Exception ex) {
        ex.printStackTrace()    
    } finally {
        sql.close()
    }

    return lastCheckedDate
}

def insertMetaRecord(String jdbcUrl, String username, String password, String taxonId) {
    def sql 
    def driver = 'com.mysql.cj.jdbc.Driver'
    def jdbcUrl = "jdbc:mysql://${params.transcriptomic_dbhost}:${params.transcriptomic_dbport}/${params.transcriptomic_dbname}"
    sql = Sql.newInstance(jdbcUrl, params.transcriptomic_dbuser,params.transcriptomic_dbpassword,driver)

    try {
        // Get the current date and time
        def currentDate = LocalDateTime.now()
        def dateFormatter = DateTimeFormatter.ofPattern("yyyy-MM-dd")
        def formattedDate = currentDate.format(dateFormatter)

        // Execute the SQL INSERT statement
        def insertQuery = "INSERT INTO meta (taxon_id, last_checked_date) VALUES ('${taxonId}', '${formattedDate}')"
        sql.executeInsert(insertQuery)
    } catch (Exception ex) {
        ex.printStackTrace()
    } finally {
        sql.close()
    }

}
def updateLastCheckedDate(String jdbcUrl, String username, String password, String taxonId) {
    def sql = Sql.newInstance(jdbcUrl, username, password)

    try {
        // Get the current date and time
        def currentDate = LocalDateTime.now()
        def dateFormatter = DateTimeFormatter.ofPattern("yyyy-MM-dd")
        def formattedDate = currentDate.format(dateFormatter)

        // Execute the SQL UPDATE statement
        def updateQuery = "UPDATE meta SET last_checked_date = '${formattedDate}' WHERE taxon_id = '${taxonId}'"
        sql.executeUpdate(updateQuery)
    } catch (Exception ex) {
        ex.printStackTrace()
    }finally {
        sql.close()
    }

}
def build_ncbi_path(gca, assembly_name) {
    final gca_splitted = gca.replaceAll("_","").tokenize(".")[0].split("(?<=\\G.{3})").join('/')
    return  'https://ftp.ncbi.nlm.nih.gov/genomes/all'  + '/' + gca_splitted + '/' + "$gca" +'_' + assembly_name.replaceAll(" ","_") + '/' + "$gca" + '_' + assembly_name.replaceAll(" ","_") + '_genomic.fna.gz'
}

def getPairedFastqsURL(String jdbcUrl, String username, String password, String run_accession) {
    def sql = Sql.newInstance(jdbcUrl, username, password)
    try {
        def query = "SELECT url FROM file INNER JOIN run ON run_id WHERE run_accession = '${run_accession}'"
        def result = sql.rows(query)
    } catch (Exception ex) {
    ex.printStackTrace()    
    } finally {
        sql.close()
    }

    return result
}

def checkFastqc(String jdbcUrl, String username, String password, String run_accession) {
    def sql = Sql.newInstance(jdbcUrl, username, password)
    def query = """ SELECT basic_statistics, per_base_sequence_quality, per_sequence_quality_scores, \
        per_base_sequence_content 
        FROM data_files df 
        INNER JOIN run r on df.run_id =r.run_id 
        WHERE r.run_id= '${run_accession}'
        """
    def qc_status = null 
    SELECT 
        df1.basic_statistics AS basic_statistics_1, df1.per_base_sequence_quality AS per_base_sequence_quality_1,
        df1.per_sequence_quality_scores AS per_sequence_quality_scores_1, df1.per_base_sequence_content AS per_base_sequence_content_1,
        df2.basic_statistics AS basic_statistics_2, df2.per_base_sequence_quality AS per_base_sequence_quality_2,
        df2.per_sequence_quality_scores AS per_sequence_quality_scores_2, df2.per_base_sequence_content AS per_base_sequence_content_2
    FROM 
        run r
    JOIN 
        data_files df1 ON r.run_id = df1.run_id
    LEFT JOIN 
        data_files df2 ON r.run_id = df2.run_id AND df1.data_file_id <> df2.data_file_id
    WHERE 
        r.run_id = '${run_accession}'
    try {
        def result = sql.rows(query)
        // Process the results
        results.each { row ->
        def basicStatistics = row.basic_statistics
        def perBaseSequenceQuality = row.per_base_sequence_quality
        def perSequenceQualityScores = row.per_sequence_quality_scores
        def perBaseSequenceContent = row.per_base_sequence_content
        if (basicStatistics=='PASS' && perBaseSequenceQuality='PASS' &&
            perSequenceQualityScores='PASS' && perBaseSequenceContent='PASS') {
            // Execute the SQL UPDATE statement
            def updateQuery = "UPDATE RUN set qc_status = 'QC_PASS' WHERE run_id= '${run_accession}'"
            sql.executeUpdate(updateQuery)
            qc_status = 'QC_PASS'
            }
        else {
            // Execute the SQL UPDATE statement
            def updateQuery = "UPDATE RUN set qc_status = 'QC_FAIL' WHERE run_id= '${run_accession}'"
            sql.executeUpdate(updateQuery)
            qc_status = 'QC_FAIL'
            }
    }
    } catch (Exception ex) {
        ex.printStackTrace()}
    finally {
        sql.close()
    }

    return qc_status
}

def checkOverrepresentedSequences(String jdbcUrl, String username, String password, String run_accession) {
    def sql = Sql.newInstance(jdbcUrl, username, password)
    def query = """ SELECT overrepresented_sequences 
        FROM data_files df 
        INNER JOIN run r on df.run_id =r.run_id 
        WHERE r.run_id= '${run_accession}'
        """
    def overrepresented_sequences = null 

    try {
        def result = sql.rows(query)
        // Process the results
        results.each { row ->
        def OverrepresentedSequences = row.overrepresented_sequences
        
        if (OverrepresentedSequences=='WARN' OR OverrepresentedSequences=='FAIL') {
            overrepresented_sequences = True
            }
        else {
            overrepresented_sequences = False
            }
    }
    } catch (Exception ex) {
        ex.printStackTrace()}
    finally {
        sql.close()
    }
    return overrepresented_sequences
}
def concatString(string1, string2, string3){
    return string1 + '_'+string2 + '_'+string3
}

def calculateIndexBases(genomeFile) {
    def indexBases = Math.min(14, Math.floor((Math.log(genomeFile, 2) / 2) - 1))
    return indexBases
}

def getRunId(String jdbcUrl, String username, String password, String run_accession, String gca, String percentage_mapped) {
    def sql = Sql.newInstance(jdbcUrl, username, password)
    def run_id = null
    try {
        def query = "SELECT run_id  FROM run WHERE run_accession = '${run_accession}'"
        run_id = sql.rows(query)
        return run_id
    } catch (Exception ex) {
        ex.printStackTrace()
    } finally {
        sql.close()
    }
    
}

def updateFastqcStatus(String jdbcUrl, String username, String password, String run_accession) {
    def sql = Sql.newInstance(jdbcUrl, username, password)

    try {
        // Execute the SQL UPDATE statement
        def updateQuery = "UPDATE run SET qc_status = 'ALIGNED' WHERE run_accession = '${run_accession}'"
        sql.executeUpdate(updateQuery)
    } catch (Exception ex) {
        ex.printStackTrace()
    }finally {
        sql.close()
    }

}