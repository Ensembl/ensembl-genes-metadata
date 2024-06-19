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
@Grab('org.codehaus.groovy:groovy-all:2.2.2')
import groovy.sql.Sql
import java.time.LocalDateTime
import java.time.format.DateTimeFormatter
import java.text.SimpleDateFormat;
import java.sql.Date

// Define global variable in the script binding
binding.driver = 'com.mysql.cj.jdbc.Driver'
binding.jdbcUrl = "jdbc:mysql://${params.transcriptomic_dbhost}:${params.transcriptomic_dbport}/${params.transcriptomic_dbname}"


def checkTaxonomy(String taxonId) {
    def sql
    sql = Sql.newInstance(jdbcUrl, params.transcriptomic_dbuser,params.transcriptomic_dbpassword,driver)
    
    try {
        def query = "SELECT * FROM meta WHERE taxon_id = ? "
        def result = sql.rows(query,[taxonId])
        return result.size() > 0
    } catch (Exception ex) {
        ex.printStackTrace()}
    finally {
        sql.close()
    }
}

def getLastCheckDate(String taxonId) {
    def sql 
    sql = Sql.newInstance(jdbcUrl, params.transcriptomic_dbuser,params.transcriptomic_dbpassword,driver)
    def result
    try {
        def query = "SELECT last_check FROM meta WHERE taxon_id = ? "
        result = sql.rows(query, [taxonId])
        if (result.size() > 0) {
            log.info("taxonomy ${taxonId} present, get last check")
            // Assuming 'last_check' is yyyy-MM-dd format
            // Adjust the date format pattern based on the actual format in your database
            //def dateFormat = new SimpleDateFormat("yyyy-MM-dd") // Adjust the format if needed
            //def dateString = Date.valueOf(result[0].last_check.toString())
            //lastCheckedDate = dateFormat.parse(dateString)

            //def dateString = result[0].last_check.toLocalDate()
            //def dateFormatter = DateTimeFormatter.ofPattern("yyyy-MM-dd")
            //lastCheckedDate = dateString.format(dateFormatter)
            lastCheckedDate = result[0].last_check
            //println(lastCheckedDate)
        }
    } catch (Exception ex) {
        ex.printStackTrace()    
    } finally {
        sql.close()
    }

    return result
}

def setLastCheckDate(String taxonId,String query_option) {
    def sql 
    sql = Sql.newInstance(jdbcUrl, params.transcriptomic_dbuser,params.transcriptomic_dbpassword,driver)

    try {
        // Get the current date and time
        def currentDate = LocalDateTime.now()
        def dateFormatter = DateTimeFormatter.ofPattern("yyyy-MM-dd")
        def formattedDate = currentDate.format(dateFormatter)
        if(query_option.contains("insert")){
        // Execute the SQL INSERT statement
        //def map = [taxon_id:'${taxonId}', last_checked_date:'${formattedDate}']
        //sql.execute "INSERT INTO meta (taxon_id, last_checked_date) VALUES ($map.taxon_id, $map.last_checked_date)"
        def params = [taxonId, formattedDate]
        sql.execute 'INSERT IGNORE INTO meta (taxon_id, last_check) VALUES (?, ?)', params
        //def insertQuery = "INSERT INTO meta (taxon_id, last_checked_date) VALUES ('${taxonId}', '${formattedDate}')"
        //sql.executeInsert(insertQuery)
        }
        if(query_option =="update"){
        // Execute the SQL INSERT statement
        def params = [formattedDate, taxonId]
        //def updateQuery = "UPDATE meta SET last_check = '${formattedDate}' WHERE taxon_id = '${taxonId}'"
        sql.execute 'UPDATE meta SET last_check = ?  WHERE taxon_id = ?', params
        }
    } catch (Exception ex) {
        ex.printStackTrace()
    } finally {
        sql.close()
    }

}

def setMetaDataRecord(String mysqlQuery){
    def sql 
    sql = Sql.newInstance(jdbcUrl, params.transcriptomic_dbuser,params.transcriptomic_dbpassword,driver)
    
    // Define a regular expression pattern to extract values within parentheses /\(([^)]+)\)/
    def pattern = /\(([^)]+)\)/
    // Match the pattern in the input string
    def matcher = (mysqlQuery =~ pattern )
    // Track the number of matches found
    def matchCount = 0
    // Iterate through matches until the second occurrence is found
    while (matcher.find()) {
        matchCount++
        if (matchCount == 2) {
        println("match found")
        // Extract the values from the match
        def paramValuesString = matcher.group(1)
        //log.info("Extracted paramValuesString: ${paramValuesString}")
        // Split the values string by commas
        println("paramValues1")
        // Trim whitespace and single quotes from each parameter value
        def paramValues = paramValuesString.split(',').collect { it.trim().replaceAll(/^'|'$/, '') }       
        // Add NULL for values that are represented as 'NULL' in the input string
        paramValues = paramValues.collect { it == 'NULL' ? null : it }.toList()
        // Print the parameter values
        log.info("Parameter values: ${paramValues}")
        // Generate a comma-separated string of placeholders
        def placeholders = paramValues.collect { "?" }.join(", ")
        // Construct the query by replacing the values string with placeholders
        def query = mysqlQuery.replaceFirst(paramValuesString, placeholders)
        // Remove the values string from the input string to get the query
        // query = mysqlQuery.replaceFirst(paramValuesString, "?".repeat(paramValues.size())).toString()
        // Print the final query
        log.info("Final query: $query")

        try{
            sql.execute query, paramValues
        } catch (Exception ex) {
            ex.printStackTrace()
        } finally {
            sql.close()
        }
    }else{
    println("no match found")
    }
}
}


def getRunTable(String runAccession, String tableKey) {
    def sql
    sql = Sql.newInstance(jdbcUrl, params.transcriptomic_dbuser,params.transcriptomic_dbpassword,driver)
    try {
        def query = "SELECT ${tableKey} FROM run WHERE run_accession = ? "
        def result = sql.firstRow(query,[runAccession])
        return result[tableKey].toString()
    } catch (Exception ex) {
        ex.printStackTrace()}
    finally {
        sql.close()
    }
}

def getDataFromTable(String queryKey, String queryTable, String tableColumn, String tableValue){
    def sql
    sql = Sql.newInstance(jdbcUrl, params.transcriptomic_dbuser,params.transcriptomic_dbpassword,driver)
    def result
    try {
        def query = "SELECT ${queryKey} FROM ${queryTable}  WHERE ${tableColumn} = ?" 
        result = sql.rows(query,[tableValue])
        
    } catch (Exception ex) {
        ex.printStackTrace()}
    finally {
        sql.close()
    }
    return result
}

def updateTable(String queryKey, String queryValue, String queryTable, String tableColumn, String tableValue) {
    def sql
    sql = Sql.newInstance(jdbcUrl, params.transcriptomic_dbuser,params.transcriptomic_dbpassword,driver)
    try {
        def query = "UPDATE ${queryTable} SET ${tableColumn}  = ? WHERE ${queryKey} = ?"
        def params = [tableValue, queryValue]
        sql.execute query, params
    } catch (Exception ex) {
        ex.printStackTrace()}
    finally {
        sql.close()
    }
}
def checkRunStatus(String runId) {
    def sql
    sql = Sql.newInstance(jdbcUrl, params.transcriptomic_dbuser,params.transcriptomic_dbpassword,driver)
    def query_fastqc_results = """
    SELECT 
        CONCAT(
            df1.basic_statistics, ', ', 
            df2.basic_statistics, ', ', 
            df1.per_base_sequence_quality, ', ',
            df2.per_base_sequence_quality, ', ', 
            df1.per_sequence_quality_scores, ', ', 
            df2.per_sequence_quality_scores, ', ',
            df1.per_base_sequence_content,  ', ', 
            df2.per_base_sequence_content
        ) AS formatted_results
    FROM 
        run r
    JOIN 
        data_files df1 ON r.run_id = df1.run_id
    JOIN 
        data_files df2 ON r.run_id = df2.run_id
    WHERE 
        r.run_id = ?
        AND df1.file_id < df2.file_id
    LIMIT 1;
    """
    def query_qc_status = """
    UPDATE run 
    SET qc_status = ?
    WHERE run_id = ? ;
    """
    try {
        def result = sql.firstRow(query_fastqc_results, [runId])
        if (result && result.formatted_results) {
            def statusList = result.formatted_results.split(', ').collect { it.trim() }
            def finalStatus = statusList.any { it == 'FAIL' } ? 'QC_FAIL' : 'QC_PASS'
            def params = [finalStatus, runId]
            sql.execute query_qc_status, params
            return finalStatus
        } else {
            return 'No data found'
        }
    } catch (Exception ex) {
        ex.printStackTrace()
        return 'Error'
    } finally {
        sql.close()
    }
}

def checkOverrepresentedSequences(String run_accession) {
    def sql = Sql.newInstance(jdbcUrl, username, password)
    def query = """ SELECT overrepresented_sequences 
                    FROM data_files df 
                    INNER JOIN run r on df.run_id =r.run_id 
                    WHERE r.run_id= ?
        """
    def overrepresented_sequences = false 

    try {
        def result = sql.rows(query,[run_accession])
        // Process the results
        results.each { row ->
        def OverrepresentedSequences = row.overrepresented_sequences
        if (OverrepresentedSequences=='WARN' || OverrepresentedSequences=='FAIL') {
            overrepresented_sequences = true
            }
        }
    } catch (Exception ex) {
        ex.printStackTrace()}
    finally {
        sql.close()
    }
    return overrepresented_sequences
}
/*
def getDataFromTable(String run_accession, String queryTable, String tableKey){
    def sql
    sql = Sql.newInstance(jdbcUrl, params.transcriptomic_dbuser,params.transcriptomic_dbpassword,driver)
    try {
        def query = "SELECT ${tableKey} FROM ${queryTable}  INNER JOIN run ON run_id WHERE run_accession = ?" 
        def result = sql.rows(query,[run_accession])
        return result.size() > 0
    } catch (Exception ex) {
        ex.printStackTrace()}
    finally {
        sql.close()
    }
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
        if (basicStatistics=='PASS' and perBaseSequenceQuality='PASS' and
            perSequenceQualityScores='PASS' and perBaseSequenceContent='PASS') {
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
*/
