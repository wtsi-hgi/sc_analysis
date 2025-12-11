#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { sc_analysis } from './workflows/sc_analysis.nf'

workflow {
    sc_analysis()
}