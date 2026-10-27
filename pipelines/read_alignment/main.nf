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
    ENSEMBL READ ALIGNMENT PIPELINE
========================================================================================

    This Nextflow DSL2 workflow performs alignment of RNA sequencing reads
    against reference genome assemblies using STAR and Minimap2.

    The pipeline supports short-read and long-read sequencing technologies,
    automatically selects the appropriate aligner based on the sequencing
    platform, and generates indexed BAM alignment files.

    Optional downstream processing includes tissue-level BAM merging,
    stranded BAM generation, BigWig coverage tracks, and CRAM conversion.

    Pipeline Stages:
    ----------------

    1. FETCH_GENOME
       Retrieve or locate reference genome assemblies specified in the
       input CSV metadata file.

    2. DOWNLOAD_FASTQS
       Download sequencing reads and prepare FASTQ files for alignment.
       Validate the presence of required FASTQ files.

    3. PLATFORM SELECTION
       Select the appropriate alignment workflow:
       - Illumina: STAR
       - PacBio / PacBio SMRT / ONT: Minimap2

    4a. STAR_INDEX_PARAMS
        Calculate genome statistics and prepare STAR indexing parameters
        for Illumina sequencing data.

    4b. STAR_INDEX_GENOME
        Generate the STAR genome index.

    4c. STAR
        Align Illumina RNA sequencing reads against the reference genome
        and produce BAM alignment files.

    5a. MINIMAP2_INDEX_GENOME
        Generate the Minimap2 reference genome index.

    5b. MINIMAP2
        Align PacBio and Oxford Nanopore reads against the reference
        genome and produce BAM alignment files.

    6. CHECK_BAM
       Validate alignment BAM files produced by STAR and Minimap2.
       Only validated BAM files proceed to downstream processing.

    7. DELETE_FASTQ
       Perform FASTQ cleanup after successful BAM validation.

    8. INDEX_BAM
       Generate BAI index files for validated BAM alignments.

    9. MERGE_BAM_PER_TISSUE
       Optionally merge BAM files belonging to the same taxon,
       sequencing platform, and tissue.

       Enabled with:
           --mergeTissue

       Merged BAM files are validated and indexed before further
       downstream processing.

    10. BAM2STRAND
        Optionally generate forward- and reverse-strand BAM files.

        Enabled with:
            --stranded

    11. BAM2BIGWIG
        Optionally generate BigWig coverage tracks.

        Enabled with:
            --bam2bigWig

        When stranded processing is enabled, BigWig generation
        receives forward- and reverse-strand BAM files.
        Otherwise, it receives the original or merged BAM files.

    12. BAM2CRAM
        Optionally convert original or merged BAM files into CRAM
        format for compressed alignment storage.

        Enabled with:
            --bam2cram

    13. INDEX_CRAM
        Generate CRAI index files for CRAM outputs.

    14. WRITE_REPORT
        Generate a final report using metadata collected from
        completed alignment and downstream processing steps.

    15. COLLECT_SOFTWARE_VERSIONS
        Collect and consolidate software version information
        from the executed modules.

    Input:
    ------

    CSV metadata file containing:
        - taxon_id
        - assembly_accession
        - platform
        - paired
        - tissue
        - run_accession
        - pair1
        - md5_1
        - pair2
        - md5_2
        - genome_file

    Rows without taxon_id, run_accession, or pair1 are filtered
    out before processing.

    The genome_file field can provide a reference genome path.
    Otherwise, the pipeline uses params.genome_file or retrieves
    the reference genome through FETCH_GENOME.

    Output:
    -------

    Depending on the enabled processing options:
        - Retrieved reference genome assemblies
        - Downloaded FASTQ files
        - STAR genome indexes
        - Minimap2 genome indexes
        - Aligned BAM files
        - BAM index files (BAI)
        - Tissue-merged BAM files
        - Forward- and reverse-strand BAM files
        - BigWig coverage tracks
        - Compressed CRAM alignment files
        - CRAM index files (CRAI)
        - Final alignment report
        - Consolidated software version information

    Execution Parameters:
    ---------------------

        --csvFile        Input CSV metadata file
        --outDir         Pipeline output directory
        --genome_file    Optional reference genome file
        --reuse_bam      Reuse existing BAM alignments when supported
        --mergeTissue    Enable tissue-level BAM merging
        --stranded       Enable stranded BAM generation
        --bam2bigWig     Enable BigWig coverage track generation
        --bam2cram       Enable BAM-to-CRAM conversion

    Workflow Notes:
    ---------------

    - STAR is used for Illumina sequencing data.
    - Minimap2 is used for PacBio, PacBio SMRT, and ONT data.
    - Both alignment branches validate and index their BAM outputs.
    - Tissue merging is optional and occurs before downstream
        stranded, BigWig, and CRAM processing.
    - CRAM conversion uses the original or tissue-merged BAM,
        independently of stranded BAM generation.
    - BigWig generation uses stranded BAMs when enabled.
    - Software versions are collected throughout the workflow.
    - SAM2BAM is imported but is not currently executed.
    - WRITE_REPORT generates the final output report.

========================================================================================
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

workflow READ_ALIGNMENT {
    take:
    csvFile

    main:
    reportInput = channel.empty()
    ch_versions_file = channel.empty()
    def finalBam
    def starOutput 
    def minimapOutput
    def bamToBigWig
    def bamForDownstream
    def data 

    data = channel.fromPath(csvFile, type: 'file', checkIfExists: true)
                .splitCsv(sep:',', header:true)
                .filter { row -> row.get('taxon_id') && row.get('run_accession') && row.get('pair1') }
                .map { row ->
                    [
                        taxon_id : row.get('taxon_id'),
                        gca : row.get('assembly_accession'),
                        platform : row.get('platform'),
                        paired : row.get('paired')?.toBoolean(),
                        tissue : row.get('tissue'),
                        run_accession : row.get('run_accession'),
                        pair1_path : row.get('pair1'),
                        md5_1 : row.get('md5_1'),
                        pair2_path : row.get('pair2') ? row.get('pair2') : null,
                        md5_2 : row.get('md5_2') ? row.get('md5_2') : null,
                        fasta_file: row.get('genome_file') ?: params.genome_file,
                        csv_path : params.csvFile,
                        output_dir : params.outDir,
                        reuse_bam : params.reuse_bam
                    ]
                }
    data.view { row -> "DATA: ${row}" }        

    def genomeAndDataToAlign = FETCH_GENOME(data).fasta_file_output
    .map {meta, genomeFile -> 
    return meta + [
        fasta_file: genomeFile,
        genome_dir: genomeFile.parent,
        alignment_dir : "${meta.output_dir}/${meta.taxon_id}/${meta.platform}/${meta.run_accession}/alignment"]}
    ch_versions_file  = ch_versions_file.mix(FETCH_GENOME.out.versions_file)

    def downloadedFastqFiles=DOWNLOAD_FASTQS(genomeAndDataToAlign).fastq_file_output
    .map { meta, reads ->
        def r1 = reads.find { it.name.endsWith('_1.fastq.gz') }
        def r2 = reads.find { it.name.endsWith('_2.fastq.gz') }

        if (!r1) error "R1 missing for ${meta.run_accession}"
        if (meta.paired && !r2) error "R2 missing for paired run ${meta.run_accession}"

        meta + [fastq1: r1, fastq2: r2]
    }
    ch_versions_file  = ch_versions_file.mix(DOWNLOAD_FASTQS.out.versions_file)
    def processReads = downloadedFastqFiles.branch { meta  ->
        def platform = meta.platform?.toString()?.toLowerCase()

        star    : platform == 'illumina'
        minimap : platform in ['pacbio', 'pacbio_smrt', 'ont']
    }
    def genomeIndexParams = STAR_INDEX_PARAMS(processReads.star).genome_stats_output.view { "PARAMS OUT: $it" }
    ch_versions_file  = ch_versions_file.mix(STAR_INDEX_PARAMS.out.versions_file)

    def genomeIndexShortData = STAR_INDEX_GENOME(genomeIndexParams).genome_index_output
    ch_versions_file  = ch_versions_file.mix(STAR_INDEX_GENOME.out.versions_file)

    def alignStarOutput = STAR(genomeIndexShortData).star_output
    .map {meta, bamFile, _spliceJunctionFile -> return tuple(meta, bamFile) }
    ch_versions_file  = ch_versions_file.mix(STAR.out.versions_file)

    def checkedBamStar = CHECK_BAM_STAR(alignStarOutput).good_bam
    ch_versions_file  = ch_versions_file.mix(CHECK_BAM_STAR.out.versions_file)

    def alignedFiles = DELETE_FASTQ_STAR(checkedBamStar).aligned_output
    ch_versions_file  = ch_versions_file.mix(DELETE_FASTQ_STAR.out.versions_file)
    
    starOutput = INDEX_BAM_STAR(alignedFiles, 'bai').indexed_output
    ch_versions_file  = ch_versions_file.mix(INDEX_BAM_STAR.out.versions_file)


    def genomeIndexLongData = MINIMAP2_INDEX_GENOME(processReads.minimap).minimap_index
    ch_versions_file  = ch_versions_file.mix(MINIMAP2_INDEX_GENOME.out.versions_file)

    def alignMinimapOutput = MINIMAP2(genomeIndexLongData).minimap_alignment
    .map {meta, bamFile -> return tuple(meta, bamFile) }
    ch_versions_file  = ch_versions_file.mix(MINIMAP2.out.versions_file)

    def checkedBamMinimap = CHECK_BAM_MINIMAP(alignMinimapOutput).good_bam
    ch_versions_file  = ch_versions_file.mix(CHECK_BAM_MINIMAP.out.versions_file)

    def cleanFile = DELETE_FASTQ_MINIMAP(checkedBamMinimap).aligned_output
    ch_versions_file  = ch_versions_file.mix(DELETE_FASTQ_MINIMAP.out.versions_file)

    //sam2bamOutput = SAM2BAM(cleanFile).sam_output
    //ch_versions_file  = ch_versions_file.mix(SAM2BAM.out.versions_file)

    minimapOutput = INDEX_BAM_MINIMAP(cleanFile, 'bai').indexed_output
    ch_versions_file  = ch_versions_file.mix(INDEX_BAM_MINIMAP.out.versions_file)            

    // Collect all aligned BAMs
    def output2process = starOutput.mix(minimapOutput)
    output2process.each { dataRow -> dataRow.view() }

    if (params.mergeTissue){
        output2process
        .map { meta, bamFile, _baiFile ->
        def groupMeta = meta.subMap('output_dir','taxon_id', 'tissue', 'platform','fasta_file','genome_dir','csv_path')
            tuple(groupMeta, bamFile)
        }
        .groupTuple()
        .view { k, v -> "GROUP: ${k} → ${v.size()} bam files" }
        .map { meta, bamFiles ->
        def newMeta = meta + [
                alignment_dir: "${meta.output_dir}/${meta.taxon_id}/${meta.platform}/${meta.tissue}/alignment"
            ]
            tuple(newMeta, bamFiles)
        }.set { bam2merge }
        bam2merge.each { dataRow -> dataRow.view() } 

        finalBam = MERGE_BAM_PER_TISSUE(bam2merge).merged_bam
        .map {meta, bamFile -> return tuple(meta, bamFile) }
        ch_versions_file  = ch_versions_file.mix(MERGE_BAM_PER_TISSUE.out.versions_file)

        checkedMergedBam = CHECK_BAM_MERGED(finalBam).good_bam
        ch_versions_file  = ch_versions_file.mix(CHECK_BAM_MERGED.out.versions_file)

        indexedMerged =  INDEX_BAM_MERGED(checkedMergedBam, 'bai')
        bamForDownstream = indexedMerged.indexed_output
        ch_versions_file  = ch_versions_file.mix(indexedMerged.versions_file)
        reportInput=reportInput.mix(indexedMerged.meta_value)
    } else{
        bamForDownstream = output2process.flatten()
        reportInput = reportInput.mix(
            bamForDownstream.map { meta, _bam -> meta })
    }

    if (params.stranded){
        bamToStrand=bamForDownstream.map{ meta, bam, _bai ->
        tuple(meta, bam)}
        strandBam = BAM2STRAND(bamToStrand)
        ch_versions_file  = ch_versions_file.mix(strandBam.versions_file)
        reportInput = reportInput.mix(strandBam.meta_value)
        bamToBigWig = strandBam.aligned_output.map { meta, forward_bam, reverse_bam, forward_index, reverse_index ->
        tuple(meta, [forward_bam, reverse_bam], [forward_index, reverse_index])
        }
    } else {
        bamToBigWig = bamForDownstream.map { meta, bam, bai ->
        tuple(meta, [bam], [bai])
        }
    }

    if(params.bam2bigWig){
        def bigWigOutput=BAM2BIGWIG(bamToBigWig)
        ch_versions_file  = ch_versions_file.mix(bigWigOutput.versions_file)
        reportInput=reportInput.mix(bigWigOutput.meta_value)
    }
    if (params.bam2cram){
        def bamToCram=bamForDownstream.map { meta, bam, _bai ->
        tuple(meta, bam)}
        def cramFile = BAM2CRAM(bamToCram)
        ch_versions_file  = ch_versions_file.mix(cramFile.versions_file)
        def indexCram = INDEX_CRAM(cramFile.cram_output,'crai')
        ch_versions_file  = ch_versions_file.mix(indexCram.versions_file)
        reportInput=reportInput.mix(indexCram.meta_value)
    }

    reportInput
        .map { meta ->
            tuple(meta.subMap('taxon_id','output_dir','csv_path'), true)
        }
        .groupTuple()
        .map { groupMeta, _doneItems ->
            groupMeta
        }
        .set { input_set }


    WRITE_REPORT(input_set)
    ch_versions_file  = ch_versions_file.mix(WRITE_REPORT.out.versions_file)

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
    READ_ALIGNMENT(params.csvFile)
}
