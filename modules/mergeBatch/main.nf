process mergeBatch {
    publishDir "${params.outdir}/",
        mode: 'copy'
    tag "QC"
    label "highMemMT1"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == "apptainer" ?
        'docker://yuanli202004/seurat5.4_doubletfinder:latest':
        'docker://yuanli202004/seurat5.4_doubletfinder:latest' }"

    input:
        tuble(batches, INFILES)

    output:
        path("${params.outdir}/QC/coreQC_all.pdf"), emit: coreQC
        path("${params.outdir}/QC/Summary_all.csv"), emit: summary

    script:
    """
    Rscript ${moduleDir}/templates/mergeBatch.R --batches ${batches} \
        --INFILES "${INFILES}" --outdir "." 
    """
}

