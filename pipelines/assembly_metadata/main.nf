#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ASSEMBLY_METADATA } from './workflows/assembly_metadata.nf'

workflow {
    ASSEMBLY_METADATA ()
}