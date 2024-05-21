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

include { setMetaDataRecord } from './utils.nf'
import java.io.File
import java.lang.StringBuffer
process STORE_METADATA {
    scratch false
    label 'python'
    tag "$run_accession"
    
    input:
    tuple val(taxon_id), val(gca), val(run_accession)
    path metadata2process
    val mysqlUpdate
    
    output:
    tuple val(taxon_id), val(gca), val(run_accession)
    stdout

    script:
    log.info("Executing Python script to get metadata for run: $run_accession")
    """

    chmod +x $projectDir/bin/write2db.py  # Set executable permissions
    write2db.py --file-path ${metadata2process} --update ${mysqlUpdate}
    """
    /*
    def process = 'python ${projectDir}/bin/write2db.py --file-path ${metadata2process} --update ${mysqlUpdate}'.execute()
    //println process.text
    process.waitFor()
    println "stdout: ${process.in.text}"
    
    def outputStream = new StringBuffer()
    def errorStream = new StringBuffer()
    process.waitForProcessOutput(outputStream, errorStream)

    // Print output and error streams
    println "Output: ${outputStream.toString()}"
    println "Error: ${errorStream.toString()}"

    
// """chmod +x $projectDir/bin/write2db.py; write2db.py --file-path ${metadata2process} --update ${mysqlUpdate} """
   // def sb = new StringBuilder()
//def proc = ['cmd','/c','echo %AAA%'].execute(["AAA=XXX", "BBB=YYY"], null)
// def proc = ['/bin/bash','-c','echo $AAA'].execute(["AAA=XXX", , "BBB=YYY"], null)
//proc.consumeProcessOutput(sb, sb)
//proc.waitForOrKill(5000)
//println sb.toString() // -> XXX #chmod +x $projectDir/bin/get_metadata.py; get_metadata.py --run ${run_accession}
//"""chmod +x $projectDir/bin/write2db.py"""
//def process = ["python3",'$projectDir/bin/write2db.py', '--file-path %file%', '--update %update%'].execute(["file=${metadata2process}", "update=${mysqlUpdate}"],null)
//def process = "python3 $projectDir/bin/write2db.py  --file-path ${metadata2process} --update ${mysqlUpdate} 2>&1".execute()
 //def cmd = ["python","$projectDir/bin/write2db.py", "--file-path", "${metadata2process}", "--update", "${mysqlUpdate}"]
 def scriptText = '''
 $projectDir/bin/write2db.py --file-path ${metadata2process}
 '''
 //]
// println cmd
    // Assign new PrintWriter to "out"
// variable of binding object.
def stringWriter = new StringWriter()
def shellBinding = new Binding(out: new PrintWriter(stringWriter))
 
// Create GroovyShell to evaluate script.
def shell = new GroovyShell(shellBinding)
 
// Run the script.
shell.evaluate(scriptText)
println stringWriter.toString()


def process = cmd.execute()
def outputStream = new StringBuffer()
def errorStream = new StringBuffer()
//def process = cmd.execute()
process.waitForProcessOutput(outputStream, errorStream)
process.waitForOrKill(1000)
    out= outputStream.toString()
    println out
    //System.err.println errorStream.toString()

    def file = new File('query.txt')
    if (!file.exists()) {
        println("File not found: query.txt")
        System.exit(1)
    }
    file.eachLine { line ->
        println(line)
    }
    
    def cmd = ["write2db.py --file-path ${metadata2process} --update ${mysqlUpdate}"]
    def process = cmd.execute()
    process.waitForOrKill(1000)
    print "Output: " + process.text
    
    log.info("Executing Python script to get metadata for run: ${metadata2process}")
    //try {
    //cmd = """python write2db.py --file-path ${metadata2process} --update ${mysqlUpdate.toString()}"""
    //def sout = new StringBuilder()
    cmd = """
    chmod +x $projectDir/bin/write2db.py; write2db.py --file-path $metadata2process --update $mysqlUpdate
    """
    proc = cmd.execute()
    // Create a buffer to capture the output
    def outputBuffer = new StringBuffer()

// Read the output of the command using an input stream
    proc.in.eachLine { line ->
    outputBuffer.append(line).append('\n')
    }
    
    // Execute the command and capture the process
    //proc = new ProcessBuilder(cmd).redirectErrorStream(true).start()

    // Read the output of the subprocess
    //output = new StringBuilder()
    //proc.inputStream.eachLine { line ->
    //    output.append(line).append('\n')
    //}
    // Wait for the subprocess to complete
    proc.waitFor()

    // Check the exit status of the subprocess
    if (proc.exitValue() != 0) {
        throw new RuntimeException("Python script failed with exit code: ${proc.exitValue()}")
    }
    output = proc.getText() // Capture the output of the process
    // Read the contents of the output file
    //def outputFile = new File('query.txt')
    //def output = outputFile.text

    // Now you can use the 'output' variable as needed
    log.info("Output of the Python script: $output")
    
    // Splitting the output into individual parts
    outputs = output.toString().split('\n')

    // Checking if there's only one output or multiple outputs
    if (outputs.size() == 1) {
        println("sono qui ")
        println(outputs)
        println(output)
        // Only one output, pass it directly
        setMetaDataRecord(outputs)
    } else {
        // Multiple outputs, loop through and pass each one
        outputs.each { singleOutput ->
            setMetaDataRecord(singleOutput)
        }
    }
    
    } catch (Exception e) {
        log.error("Error executing Python script: ${e.message}")
        throw e // Rethrow the exception to halt the process
    }
   */ 
}
