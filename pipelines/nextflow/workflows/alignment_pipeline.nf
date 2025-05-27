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


include { FETCH_GENOME } from '../modules/fetch_genome'
include { DOWNLOAD_FASTQS } from '../modules/download_fastqs'
include { STAR_INDEX_GENOME } from '../modules/star_index_genome'
include { STAR } from '../modules/star'
include { MINIMAP2_INDEX_GENOME } from '../modules/minimap2_index_genome'
include { MINIMAP2 } from '../modules/minimap2'
include { SAM2BAM } from '../modules/sam2bam'
include { MERGE_BAM_PER_TISSUE } from '../modules/merge_bam_per_tissue'
include { BAM2STRAND } from '../modules/bam2strand'
include { INDEXING_FILES as INDEX_BAM } from '../modules/indexing_files'
include { INDEXING_FILES as INDEX_CRAM } from '../modules/indexing_files'
include { INDEXING_FILES as INDEX_BIGWIG } from '../modules/indexing_files'
include { BAM2CRAM } from '../modules/bam2cram'
include { BAM2BIGWIG } from '../modules/bam2bigWig'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ALIGNMENT_PIPELINE {
    take:
    csvFile
    bam2cram
    mergeTissue
    stranded
    bam2bigWig

    main:
    def data = Channel.fromList(samplesheetToList(csvFile, "${projectDir}/assets/schema_input.json"))
            .map { row -> 
                def url2 = row.paired == "true" ? row.pair2 : null
                def md5_2 = row.paired == "true" ? row.md5_2 : null
                tuple(taxon_id:row.get('taxon_id'), gca:row.get('gca'),platform:row.get('platform'), \
                paired:row.get('paired'),
                tissue:row.get('tissue'),run_accession:row.get('run_accession'), pair1:row.get('pair1'),
                md5_1:row.get('md5_1'),pair2:url2,md5_2:md5_2)}
    //def data = Channel.fromPath(params.csvFile, type: 'file', checkIfExists: true)
    //           .splitCsv(sep:',', header:true)
    //            .map { row -> [taxon_id:row.get('taxon_id'), gca:row.get('gca'), \
    //            run_accession:row.get('run_accession'), pair1:row.get('pair1'),
    //            pair2:row.get('pair2'),md5_1:row.get('md5_1'),md5_2:row.get('md5_2')]}
    data.each { dataRow -> dataRow.view() }    
    def alignOutput
    def finalBam

    def genomeAndDataToAlign = FETCH_GENOME(data.flatten())
    def downloadedFastqFiles=DOWNLOAD_FASTQS(genomeAndDataToAlign)
    if (downloadedFastqFiles.paired=='true'){
        def genomeIndexShortData = STAR_INDEX_GENOME(downloadedFastqFiles)
        alignStarOutput = STAR(genomeIndexShortData)
        alignOutput = INDEX_BAM (alignStarOutput,'bai')
    } else {
        def genomeIndexLongData=MINIMAP2_INDEX_GENOME(downloadedFastqFiles)
        def minimapOutput = MINIMAP2(genomeIndexLongData)
        alignOutput = SAM2BAM(minimapOutput)
    }
    // Collect all aligned BAMs
    def output2process = alignOutput.collect()
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
    if (mergeTissue){
        // Channel of [taxon_id, genomeDir, tissue, run_accession, bamFile]
        Channel.from(output2process)
        // Group by the 3rd element: tissue
        .groupBy { _taxon_id, _genomeDir, tissue, _run_accession, _bamFile -> tissue }
        // For each tissue, pull out taxon_id, genomeDir from the first entry, drop run_accession,
        // and collect all the bamFile paths into a List<String>
        .map { tissue, samples ->
            // entries is a List of 5-element lists
            def (firstTaxonId, firstGenomeDir, _el2, _el3, _el4) = samples[0]
            def bamFiles = samples.collect { it[4] }
            // emit exactly 4 elements: taxon_id, genomeDir, tissue, bamFiles
            tuple(
                firstTaxonId,
                firstGenomeDir,
                tissue,
                tissue,
                bamFiles
            )
        }
        .set { bam2merge }

        finalBam = MERGE_BAM_PER_TISSUE(bam2merge).flatten()
        //finalBam = output2process.groupBy{ it.tissue }.flatten()

    } else{
        //finalBam = output2process.flatten()
        Channel.fromList(output2process).set { finalBam }
    }
    if (stranded){
        def bamToStrand=finalBam
        BAM2STRAND(bamToStrand) //strandedBamOutput = 
        //def indexBamOutput = INDEX_BAM (strandedBamOutput,'bai')
    }
    if (bam2cram){
        def bamToCram=finalBam
        def cramFile = BAM2CRAM (bamToCram)
        INDEX_CRAM (cramFile,'crai') //indexCramFile
    }
    if(bam2bigWig){
        def bamToBigWig=finalBam
        BAM2BIGWIG(bamToBigWig) //bigWigOutput

    }


}  


