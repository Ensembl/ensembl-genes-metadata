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

nextflow.enable.dsl=2



// Load Plugins
include { validateParameters; paramsSummaryLog; samplesheetToList} from 'plugin/nf-schema'




/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { FETCH_GENOME } from './modules/fetch_genome.nf'
include { DOWNLOAD_FASTQS } from './modules/download_fastqs.nf'
include { STAR_INDEX_GENOME } from './modules/star_index_genome.nf'
include { STAR } from './modules/star.nf'
include { MINIMAP2_INDEX_GENOME } from './modules/minimap2_index_genome.nf'
include { MINIMAP2 } from './modules/minimap2.nf'
include { SAM2BAM } from './modules/sam2bam.nf'
include { MERGE_BAM_PER_TISSUE } from './modules/merge_bam_per_tissue.nf'
include { BAM2STRAND } from './modules/bam2strand.nf'
include { INDEXING_FILES as INDEX_BAM } from './modules/indexing_files.nf'
include { INDEXING_FILES as INDEX_CRAM } from './modules/indexing_files.nf'
include { INDEXING_FILES as INDEX_BIGWIG } from './modules/indexing_files.nf'
include { BAM2CRAM } from './modules/bam2cram.nf'
include { BAM2BIGWIG } from './modules/bam2bigWig.nf'
include { DELETE_FASTQ } from './modules/delete_fastq.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    take:
    csvFile
    bam2cram
    mergeTissue
    stranded
    bam2bigWig

    main:
    // Validate input parameters
    validateParameters()

    // Print summary of supplied parameters
    log.info(paramsSummaryLog(workflow))
    ch_versions_file = channel.empty()
    paired_sample = ""
    def data 
    data = channel.fromPath(csvFile, type: 'file', checkIfExists: true)
                .splitCsv(sep:',', header:true)
                .filter { row -> row.get('taxon_id') && row.get('run_accession') && row.get('pair1') }
                .map { row ->
                    [
                        taxonId : row.get('taxon_id'),
                        assembly_accession : row.get('gca'),
                        instrument_platform : row.get('platform'),
                        paired : row.get('paired')?.toBoolean(),
                        tissue_name : row.get('tissue'),
                        run_accession : row.get('run_accession'),
                        pair1_path : row.get('pair1'),
                        md5_1_sum : row.get('md5_1'),
                        pair2_path : paired_sample ? row.get('pair2') : null,
                        md5_2 : paired_sample ? row.get('md5_2') : null,
                        genome_file: row.get('genome_file')
                    ]
                }
    data.view { row -> "DATA: ${row}" }        

    def alignOutput
    def finalBam

    def genomeAndDataToAlign = FETCH_GENOME(data).fasta_file_output
    .map {meta, genomeFile -> return tuple(meta + [fasta_file: genomeFile]) }
    ch_versions_file = ch_versions_file.mix(FETCH_GENOME.out.versions_file)

    def downloadedFastqFiles=DOWNLOAD_FASTQS(genomeAndDataToAlign).fastq_file_output
    .map {meta, fastq1, fastq2 -> return tuple(meta + [fastq1: fastq1, fastq2: fastq2]) }
    ch_versions_file = ch_versions_file.mix(DOWNLOAD_FASTQS.out.versions_file)
    if(downloadedFastqFiles.filter { it[4] == true }){
        log.info("Paired-end data detected. Using STAR for alignment.")
    //|| downloadedFastqFiles[5].toString() == 'true') {
    //if (downloadedFastqFiles[4]==true)
        def genomeIndexShortData = STAR_INDEX_GENOME(downloadedFastqFiles)
        alignStarOutput = STAR(genomeIndexShortData)
        deleteZipfiles = DELETE_FASTQ(alignStarOutput)
        alignOutput = INDEX_BAM (deleteZipfiles,'bai')
        //.map { result ->
     //def (taxon_id, genomeDir, tissue, platform, run_accession, bamFile) = result
     //def key = [taxon_id, genomeDir, platform, tissue]
      //      tuple(key, [run_accession: run_accession, bamFile: bamFile])
    //}.collect()

        //.map { result ->
   // def (taxon_id, genomeDir, tissue, platform, run_accession, bamFile) = result
   //return [taxon_id, genomeDir, tissue, platform, run_accession, bamFile]} 
   alignOutput.each { dataRow -> dataRow.view() }
   } else {
   downloadedFastqFiles.each { dataRow -> dataRow[5].view() }    
        def genomeIndexLongData=MINIMAP2_INDEX_GENOME(downloadedFastqFiles)
        def minimapOutput = MINIMAP2(genomeIndexLongData)
        cleanFile = DELETE_FASTQ(minimapOutput)
        alignOutput = SAM2BAM(cleanFile)
    }
    // Collect all aligned BAMs
    def output2process = alignOutput
    //.map { result ->
    //def (taxon_id, genomeDir, gca, platform, paired, tissue, run_accession, bamFile) = result
    //return [
    //    taxon_id    : taxon_id,
    //    genomeDir    : genomeDir,
    //    gca         : gca,
    //    platform    : platform,
    //    paired      : paired,
    //    tissue      : tissue,
    //    run_accession: run_accession,
    //   bamFile         : bamFile
    //]
    //    }
    
    //output_dir in the following processes will be tissue if merging is eanbled otherwise it will be the run_accession
    //def output2process = alignOutput
//    def output2process = alignOutput.collect().map { result ->
  //   def (taxon_id, genomeDir, tissue, platform, run_accession, bamFile) = result
   // [taxon_id, genomeDir, tissue, platform, run_accession, bamFile]
//    def key = [taxon_id, genomeDir, platform, tissue]
  //          tuple(key, [run_accession: run_accession, bamFile: bamFile])
    //}
    //log.info("data: ${data}")
    //return [taxon_id:taxon_id, genomeDir:genomeDir, tissue:tissue, platform:platform, run_accession:run_accession, bamFile:bamFile]
    //}.collect()
//    return [taxon_id, genomeDir, tissue, platform, run_accession, bamFile]
//}.toList()
    output2process.each { dataRow -> dataRow.view() }
    if (mergeTissue){
        // Group by the 3rd element: tissue
        //.groupBy { _taxon_id, _genomeDir, tissue, _run_accession, _bamFile -> tissue }
        //.map { tissue, samples ->
            // entries is a List of 5-element lists
        //    def (firstTaxonId, firstGenomeDir, _el2, _el3, _el4) = samples[0]
         //   def bamFiles = samples.collect { it[4] }
            // emit exactly 4 elements: taxon_id, genomeDir, tissue, bamFiles
           // tuple(
           //     firstTaxonId,
          //      firstGenomeDir,
          //      tissue,
           //     tissue,
           //     bamFiles
           // )
       // }
       // .set { bam2merge }
//       FULL_LENGTH_CONTAINS.out.full_length_contains
 //      .map { _chromosome_id, genome_id, _full_length_contains, full_length_contains_stats, chromosome_count -> 
   //    tuple(groupKey(genome_id, chromosome_count.toInteger())), full_length_contains_stats }.groupTuple()
   output2process
    .map { row ->
        def (taxon_id, genomeDir, tissue, platform, output_dir, bamFile) = row
        return [tuple(taxon_id, genomeDir, tissue, platform),bamFile]
    }
    .groupTuple()  // group by (taxon_id, tissue)
   .view { k, v -> "GROUP: ${k} → ${v.size()} bam files" }
   .map { groupKey, values ->
        def (taxon_id, genomeDir, tissue, platform) = groupKey
        //def genomeDir = values[0][0]  // assume same genomeDir for group
        def bamFiles  = values.collect { it }  // collect all bamFiles
        tuple(taxon_id, genomeDir, tissue, platform, bamFiles)
    }
    .set { bam2merge }
bam2merge.each { dataRow -> dataRow.view() } 

        finalBam = MERGE_BAM_PER_TISSUE(bam2merge)
        //finalBam = output2process.groupBy{ it.tissue }.flatten()
finalBam.each { dataRow -> dataRow.view() }
    } else{
        finalBam = output2process.flatten()
    
        //Channel.fromList(output2process).set { finalBam }
    }
    if (stranded){
        def bamToStrand=finalBam
        def strandOutput=BAM2STRAND(bamToStrand) //strandedBamOutput = 
        //def indexBamOutput = INDEX_BAM (strandedBamOutput,'bai')
        if(bam2bigWig){
        BAM2BIGWIG(strandOutput) //bigWigOutput
        }

    } else {
if(bam2bigWig){
        //def bamToBigWig=finalBam
        finalBam.map { row ->
        def (taxon_id, genomeDir, tissue, platform, output_dir, bamFile) = row
        return [taxon_id, genomeDir, tissue, platform,output_dir,bamFile,file("dummy.bam")]
    }.set { bamToBigWig }
        BAM2BIGWIG(bamToBigWig) //bigWigOutput

    }
    }
    if (bam2cram){
        def bamToCram=finalBam
        def cramFile = BAM2CRAM (bamToCram).map {row -> 
        def(taxon_id, genomeDir, tissue, platform, cram_file) = row 
        def output_dir="${platform}/${tissue}"
        return tuple(taxon_id, genomeDir, tissue, platform, output_dir, cram_file)
    }//val(taxon_id), val(genomeDir), val(tissue), val(platform), path("*.cram")
        INDEX_CRAM (cramFile,'crai') //indexCramFile
    }
    //if(bam2bigWig){
    //    def bamToBigWig=finalBam
    //    BAM2BIGWIG(bamToBigWig) //bigWigOutput

    //}
    if( !stranded && !bam2cram && !bam2bigWig ) {
    println "❌ No processing options selected (stranded, bam2cram, bam2bigWig)."
}


}  


