#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cosMxDS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : 
    Contact: Yuan.li@nbis.se
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

println """\
         cosMxDS   P I P E L I N E
         ===================================
         GitHub: 
         ___________________________________
         OUTPUT DIR     : ${params.outdir}
        ___________________________________
         """
         .stripIndent()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include { FILTERBAM } from '../modules/FILTERBAM'
include { QC_workflow } from './subworkflows/1QC.nf'
include { mergeBatch } from './modules/3mergeBatch'

workflow ENTRY_QC {
    batch = Channel.from(params.batch.split(','))
    QC_workflow(batch)
}

workflow ENTRY_mergeData {
    batches = Channel.value(params.batch)
    INFILES = Channel.value(params.INFILES)
    ch_combined = batches.combine(INFILES)
    mergeBatch(ch_combined)
}

workflow {
    // QC
    batch = Channel.from(params.batch.split(','))
    QC_workflow(batch)

}

workflow.onComplete {
    println( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
    )
}

