#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ASSEMBLY_METADATA_UPDATE } from './workflow/assembly_metadata_update.nf'

workflow {
    ASSEMBLY_METADATA_UPDATE ()
}