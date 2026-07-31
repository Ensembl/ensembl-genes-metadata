#!/usr/bin/env nextflow

include { validateParameters ; paramsSummaryLog } from 'plugin/nf-schema'
include { ASSEMBLY_METADATA_UPDATE } from './workflow/assembly_metadata_update.nf'

workflow {
    validateParameters()
    log.info(paramsSummaryLog(workflow))

    if (params.gca_input && params.screen_date) {
        error("--screen_date and --gca_input are mutually exclusive. Use --gca_input with --gca_list to check a specific list of accessions, or --screen_date to screen the database.")
    }
    if (!params.gca_input && !params.screen_date) {
        error("Either --screen_date or --gca_input (with --gca_list) must be provided.")
    }

    ASSEMBLY_METADATA_UPDATE()
}
