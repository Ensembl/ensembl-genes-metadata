#!/usr/bin/env nextflow

include { validateParameters ; paramsSummaryLog } from 'plugin/nf-schema'
include { ASSEMBLY_METADATA_UPDATE } from './workflow/assembly_metadata_update.nf'

workflow {
    validateParameters()
    log.info(paramsSummaryLog(workflow))

    ASSEMBLY_METADATA_UPDATE()
}
