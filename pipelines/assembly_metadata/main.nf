#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
include { ASSEMBLY_METADATA } from './workflow/assembly_metadata.nf'

workflow {
    validateParameters()
    log.info(paramsSummaryLog(workflow))

    ASSEMBLY_METADATA ()
}