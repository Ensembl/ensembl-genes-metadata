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
import java.util.concurrent.TimeUnit

// Define global variable in the script binding
binding.driver = 'com.mysql.cj.jdbc.Driver'
binding.jdbcUrl = "jdbc:mysql://${params.transcriptomic_dbhost}:${params.transcriptomic_dbport}/${params.transcriptomic_dbname}"

class RetryUtil {
    static def retry(int maxAttempts, long delayMillis, Closure closure) {
        int attempt = 0
        while (attempt < maxAttempts) {
            try {
                return closure.call()
            } catch (Exception e) {
                attempt++
                if (attempt >= maxAttempts) {
                    throw e
                }
                TimeUnit.MILLISECONDS.sleep(delayMillis)
            }
        }
    }
}
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
            lastCheckedDate = result[0].last_check
            log.info("taxonomy ${taxonId} present, last check on ${lastCheckedDate}")
            // Assuming 'last_check' is yyyy-MM-dd format
            // Adjust the date format pattern based on the actual format in your database
            //def dateFormat = new SimpleDateFormat("yyyy-MM-dd") // Adjust the format if needed
            //def dateString = Date.valueOf(result[0].last_check.toString())
            //lastCheckedDate = dateFormat.parse(dateString)
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
        def firstDate = Date.valueOf("2019-01-01")
        if(query_option.contains("insert")){
        // Execute the SQL INSERT statement
        def params = [taxonId, firstDate]
        sql.execute 'INSERT IGNORE INTO meta (taxon_id, last_check) VALUES (?, ?)', params
        }
        if(query_option =="update"){
        // Execute the SQL INSERT statement
        def params = [formattedDate, taxonId]
        sql.execute 'UPDATE meta SET last_check = ?  WHERE taxon_id = ?', params
        }
    } catch (Exception ex) {
        ex.printStackTrace()
    } finally {
        sql.close()
    }

}

def setMetaDataRecord(String mysqlQuery){
  //  synchronized(this) {
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
        // Extract the values from the match
        def paramValuesString = matcher.group(1)
        //log.info("Extracted paramValuesString: ${paramValuesString}")
        // Split the values string by commas
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
        // Print the final query
        log.info("Final query: $query")

        try{
            RetryUtil.retry(5, 1000) {   
                    sql.withTransaction {
                        sql.execute query, paramValues
                }
            } 
        } catch (Exception ex) {
            ex.printStackTrace()
        } finally {
            sql.close()
        }
    }else{
    log.info("no match found")
    }
}
}


def getDataFromTable(String queryKey, String queryTable, String tableColumn, String tableValue){
    //def sql
    //sql = Sql.newInstance(jdbcUrl, params.transcriptomic_dbuser,params.transcriptomic_dbpassword,driver)
    def result
    def retries = 5 // Set the number of retries
    def retryCount = 0
    def sleepTime = 2000 // Set the initial sleep time in milliseconds (2 seconds)

    while (retryCount < retries) {
        def sql
        try {
            sql = Sql.newInstance(jdbcUrl, params.transcriptomic_dbuser,params.transcriptomic_dbpassword,driver)
            def query = "SELECT ${queryKey} FROM ${queryTable}  WHERE ${tableColumn} = ?" 
            result = sql.rows(query,[tableValue])
            if (result && !result.isEmpty()) {
                return result // Return the result if it's not null
            } else {
                println("No result found, retrying... (${retryCount + 1}/${retries})")
                retryCount++
                Thread.sleep(sleepTime) // Wait before retrying
            }                                                                                            
        
    } catch (Exception ex) {
        ex.printStackTrace()
        retryCount++
        Thread.sleep(sleepTime) // Wait before retrying
    }finally {
        sql.close()
    }
    }
    //return result
    println("Failed to retrieve data after ${retries} attempts")
    return 
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
def checkRunFastQCStatus(String runId) {
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
            log.info("statusList ${statusList}")
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

def calculateIndexBases(genomeFile) {
    def indexBases = Math.min(14, Math.floor((Math.log(genomeFile, 2) / 2) - 1))
    return indexBases
}
*/
def deleteRecursively(Path path) {
    if (Files.isDirectory(path)) {
        Files.newDirectoryStream(path).each { subPath ->
            deleteRecursively(subPath)
        }
    }
    Files.delete(path)
    println "Deleted: ${path}"
}
