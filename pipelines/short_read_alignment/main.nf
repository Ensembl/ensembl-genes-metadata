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

/*
========================================================================================
    SHORT READ ALIGNMENT PIPELINE
========================================================================================
    This workflow performs short read alignment on genome assemblies using
    STAR and Minimap2 for alignment.
    

    Pipeline Stages:
    1. FETCH_GENOME           - Download genome assemblies from NCBI
    2. DOWNLOAD_FASTQS        - Download FASTQ files for alignment
    3. STAR_INDEX_GENOME      - Build STAR index for genome
    4. MINIMAP2_INDEX_GENOME  - Build Minimap2 index for genome
    5. STAR                   - Perform alignment using STAR
    6. MINIMAP2               - Perform alignment using Minimap2
    7. SAM2BAM                 - Convert SAM files to BAM format
    8. MERGE_BAM_PER_TISSUE   - Merge BAM files per tissue
    9. BAM2STRAND               - Convert BAM files to stranded format
    10. INDEXING_FILES          - Index BAM, CRAM, and BigWig
    11. BAM2CRAM                 - Convert BAM files to CRAM format
    12. BAM2BIGWIG               - Convert BAM files to BigWig format
    13. DELETE_FASTQ              - Delete FASTQ files after processing
    14. COLLECT_SOFTWARE_VERSIONS - Collect software versions from all processes and merge them into a single file
    Input:
        CSV file see schema for details
    
    Output:
        -short read alignment files in BAM, CRAM, and BigWig formats
        -final report with the path for all output files
----------------------------------------------------------------------------------------
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
include { STAR_INDEX_PARAMS } from './modules/star_index_params.nf'
include { STAR } from './modules/star.nf'
include { MINIMAP2_INDEX_GENOME } from './modules/minimap2_index_genome.nf'
include { MINIMAP2 } from './modules/minimap2.nf'
include { SAM2BAM } from './modules/sam2bam.nf'
include { MERGE_BAM_PER_TISSUE } from './modules/merge_bam_per_tissue.nf'
include { BAM2STRAND } from './modules/bam2strand.nf'
include { INDEXING_FILES as INDEX_CRAM } from './modules/indexing_files.nf'
include { INDEXING_FILES as INDEX_BIGWIG } from './modules/indexing_files.nf'
include { INDEXING_FILES as INDEX_BAM_STAR } from './modules/indexing_files.nf'
include { INDEXING_FILES as INDEX_BAM_MINIMAP } from './modules/indexing_files.nf'
include { INDEXING_FILES as INDEX_BAM_MERGED } from './modules/indexing_files.nf'
include { BAM2CRAM } from './modules/bam2cram.nf'
include { BAM2BIGWIG } from './modules/bam2bigWig.nf'
include { DELETE_FASTQ as DELETE_FASTQ_STAR } from './modules/delete_fastq.nf'
include { DELETE_FASTQ as DELETE_FASTQ_MINIMAP } from './modules/delete_fastq.nf'
include { CHECK_BAM as CHECK_BAM_STAR } from './modules/check_bam.nf'
include { CHECK_BAM as CHECK_BAM_MINIMAP } from './modules/check_bam.nf'
include { CHECK_BAM as CHECK_BAM_MERGED } from './modules/check_bam.nf'
include { COLLECT_SOFTWARE_VERSIONS } from './modules/collect_software_versions.nf'
include { WRITE_REPORT } from './modules/write_report.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SHORT_READ_ALIGNMENT {
    take:
    csvFile

    main:
    ch_versions_file = channel.empty()
    paired_sample = ""
    def data 
    data = channel.fromPath(csvFile, type: 'file', checkIfExists: true)
                .splitCsv(sep:',', header:true)
                .filter { row -> row.get('taxon_id') && row.get('run_accession') && row.get('pair1') }
                .map { row ->
                    [
                        taxonId : row.get('taxon_id'),
                        gca : row.get('assembly_accession'),
                        instrument_platform : row.get('platform'),
                        paired : row.get('paired')?.toBoolean(),
                        tissue_name : row.get('tissue'),
                        run_accession : row.get('run_accession'),
                        pair1_path : row.get('pair1'),
                        md5_1_sum : row.get('md5_1'),
                        pair2_path : paired_sample ? row.get('pair2') : null,
                        md5_2 : paired_sample ? row.get('md5_2') : null,
                        fasta_file: row.get('genome_file') || params.genome_file
                    ]
                }
    data.view { row -> "DATA: ${row}" }        

    //def alignOutput
    def finalBam
    def starOutput 
    def minimapOutput

    def genomeAndDataToAlign = FETCH_GENOME(data).fasta_file_output
    .map {meta, genomeFile -> return [meta + [fasta_file: genomeFile]] }
    ch_versions_file = ch_versions_file.mix(FETCH_GENOME.out.versions_file)

    def downloadedFastqFiles=DOWNLOAD_FASTQS(genomeAndDataToAlign).fastq_file_output
    .map {meta, fastq1, fastq2 -> return [meta + [fastq1: fastq1, fastq2: fastq2]] }
    ch_versions_file = ch_versions_file.mix(DOWNLOAD_FASTQS.out.versions_file)
    def processReads = downloadedFastqFiles.branch { meta  ->
        def platform = meta.instrument_platform?.toString()?.toLowerCase()

        star    : platform == 'illumina'
        minimap : platform in ['pacbio', 'pacbio_smrt', 'ont']
    }
    //if (processReads.star) {
        //log.info("Illumina data detected. Using STAR for alignment.")
        def genomeIndexParams = STAR_INDEX_PARAMS(processReads.star).genome_stats_output
        ch_versions_file = ch_versions_file.mix(STAR_INDEX_PARAMS.out.versions_file)

        def genomeIndexShortData = STAR_INDEX_GENOME(genomeIndexParams).genome_index_output
        .map {meta, genomeDir -> return [meta + [genome_dir: genomeDir]] }
        ch_versions_file = ch_versions_file.mix(STAR_INDEX_GENOME.out.versions_file)

        def alignStarOutput = STAR(genomeIndexShortData).star_output
        .map {meta, output_dir, bamFile -> return tuple(meta + [output_dir: output_dir], bamFile) }
        ch_versions_file = ch_versions_file.mix(STAR.out.versions_file)

        def checkedBamStar = CHECK_BAM_STAR(alignStarOutput).good_bam
        ch_versions_file = ch_versions_file.mix(CHECK_BAM_STAR.out.versions_file)

        def alignedFiles = DELETE_FASTQ_STAR(checkedBamStar).aligned_output
        ch_versions_file = ch_versions_file.mix(DELETE_FASTQ_STAR.out.versions_file)
        
        starOutput = INDEX_BAM_STAR(alignedFiles, 'bai').aligned_output
        ch_versions_file = ch_versions_file.mix(INDEX_BAM_STAR.out.versions_file)
    //}

    //if (processReads.minimap) {
        //log.info("PacBio data detected. Using Minimap2 for alignment.")

        def genomeIndexLongData = MINIMAP2_INDEX_GENOME(processReads.minimap).minimap_index
        ch_versions_file = ch_versions_file.mix(MINIMAP2_INDEX_GENOME.out.versions_file)

        def alignMinimapOutput = MINIMAP2(genomeIndexLongData).minimap_alignment
        .map {meta, output_dir, bamFile -> return tuple(meta + [output_dir: output_dir], bamFile) }
        ch_versions_file = ch_versions_file.mix(MINIMAP2.out.versions_file)

        def checkedBamMinimap = CHECK_BAM_MINIMAP(alignMinimapOutput).good_bam
        ch_versions_file = ch_versions_file.mix(CHECK_BAM_MINIMAP.out.versions_file)

        def cleanFile = DELETE_FASTQ_MINIMAP(checkedBamMinimap).aligned_output
        ch_versions_file = ch_versions_file.mix(DELETE_FASTQ_MINIMAP.out.versions_file)

        sam2bamOutput = SAM2BAM(cleanFile).sam_output
        ch_versions_file = ch_versions_file.mix(SAM2BAM.out.versions_file)

        minimapOutput = INDEX_BAM_MINIMAP(sam2bamOutput, 'bai').aligned_output
        ch_versions_file = ch_versions_file.mix(INDEX_BAM_MINIMAP.out.versions_file)
    //}            

    // Collect all aligned BAMs
    def output2process = starOutput.mix(minimapOutput)
    output2process.each { dataRow -> dataRow.view() }
    def mergedBam
    if (params.mergeTissue){
        output2process
        .map { meta, bamFile ->
            def groupMeta = meta.subMap('taxonId', 'tissue_name', 'platform')
            tuple(groupMeta, bamFile)
        }
        .groupTuple()
        .view { k, v -> "GROUP: ${k} → ${v.size()} bam files" }
        .map { meta, bamFiles ->
            tuple(meta, bamFiles)
        }.set { bam2merge }
    //output2process.each { dataRow -> dataRow.view() }
    //if (mergeTissue){
    //output2process
    //.map { row ->
    //    def (taxon_id, genomeDir, tissue, platform, output_dir, bamFile) = row
    //    return [tuple(taxon_id, genomeDir, tissue, platform),bamFile]
    //}
    //.groupTuple()  // group by (taxon_id, tissue)
    //.view { k, v -> "GROUP: ${k} → ${v.size()} bam files" }
    //.map { groupKey, values ->
    //        def (taxon_id, genomeDir, tissue, platform) = groupKey
            //def genomeDir = values[0][0]  // assume same genomeDir for group
    //        def bamFiles  = values.collect { it }  // collect all bamFiles
    //        tuple(taxon_id, genomeDir, tissue, platform, bamFiles)
    //    }
    //        .set { bam2merge }
        bam2merge.each { dataRow -> dataRow.view() } 

        finalBam = MERGE_BAM_PER_TISSUE(bam2merge).merged_bam
        .map {meta, output_dir, bamFile -> return tuple(meta + [output_dir: output_dir], bamFile) }
        ch_versions_file = ch_versions_file.mix(MERGE_BAM_PER_TISSUE.out.versions_file)

        checkedMergedBam = CHECK_BAM_MERGED(finalBam).good_bam
        ch_versions_file = ch_versions_file.mix(CHECK_BAM_MERGED.out.versions_file)

        mergedBam = INDEX_BAM_MERGED(checkedMergedBam, 'bai').aligned_output
        ch_versions_file = ch_versions_file.mix(INDEX_BAM_MERGED.out.versions_file)
        mergedBam.each { dataRow -> dataRow.view() }

        } else{
            finalBam = output2process.flatten()
        }
    //Define a finalBam channel to hold the final BAM files after merging or flattening
    def bamForDownstream = params.mergeTissue ? mergedBam : output2process        
    if (params.stranded){
        def bamToStrand=bamForDownstream
        def strandOutput=BAM2STRAND(bamToStrand).aligned_output
        ch_versions_file = ch_versions_file.mix(BAM2STRAND.out.versions_file)

        if(params.bam2bigWig){
            BAM2BIGWIG(strandOutput)
            ch_versions_file = ch_versions_file.mix(BAM2BIGWIG.out.versions_file)
        }
    
    } else {
    if(params.bam2bigWig){
            def bamToBigWig=bamForDownstream
            //bamForDownstream.map { row ->
            //def (taxon_id, genomeDir, tissue, platform, output_dir, bamFile) = row
            //return [taxon_id, genomeDir, tissue, platform,output_dir,bamFile,file("dummy.bam")]
        //}.set { bamToBigWig }
            BAM2BIGWIG(bamToBigWig)
            ch_versions_file = ch_versions_file.mix(BAM2BIGWIG.out.versions_file)
        }
    }
    if (params.bam2cram){
        def bamToCram=bamForDownstream
        def cramFile = BAM2CRAM(bamToCram).cram_output
        ch_versions_file = ch_versions_file.mix(BAM2CRAM.out.versions_file)
        //.map {row -> 
        //def(taxon_id, genomeDir, tissue, platform, cram_file) = row 
        //def output_dir="${platform}/${tissue}"
        //return tuple(taxon_id, genomeDir, tissue, platform, output_dir, cram_file)
        //}
        INDEX_CRAM (cramFile,'crai') //indexCramFile
        ch_versions_file = ch_versions_file.mix(INDEX_CRAM.out.versions_file)
    }
    def reportInput = bamForDownstream
    .map { meta, bam -> meta }
    .collect()
    WRITE_REPORT(reportInput)
    ch_versions_file = ch_versions_file.mix(WRITE_REPORT.out.versions_file)

    // Merge into single file and publish
    COLLECT_SOFTWARE_VERSIONS(ch_versions_file.collect())
    if( !params.stranded && !params.bam2cram && !params.bam2bigWig ) {
    println "❌ No processing options selected (stranded, bam2cram, bam2bigWig)."
}


}  

workflow {
    log.info("Pipeline started at: ${new Date().format('dd-MM-yyyy HH:mm:ss')}")
    // Validate input parameters
    validateParameters()
    // Print summary of supplied parameters
    log.info(paramsSummaryLog(workflow))
    // Execute main workflow
    SHORT_READ_ALIGNMENT(params.csvFile)
}
