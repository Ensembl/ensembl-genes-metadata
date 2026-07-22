#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
include { ASSEMBLY_METADATA_UPDATE } from './workflow/assembly_metadata_update.nf'

validateParameters()

log.info paramsSummaryLog(workflow)

workflow {
    ASSEMBLY_METADATA_UPDATE ()
}