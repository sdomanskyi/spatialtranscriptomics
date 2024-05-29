#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ST } from './workflows/spatialtranscriptomics'

workflow {

    Channel
    .from(file(params.input))
    .splitCsv(header:true, sep:',')
    .map( { it -> [ (it.sample_id), it ] } )
    .set{ samples }

    ST ( samples )

}
