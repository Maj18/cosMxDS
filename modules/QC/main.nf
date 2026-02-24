process QC {
    publishDir "${params.outdir}/",
        mode: 'copy'
    tag "QC"
    label "highMemMT2"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == "apptainer" ?
        'docker://yuanli202004/seurat5.4_doubletfinder:latest':
        'docker://yuanli202004/seurat5.4_doubletfinder:latest' }"

    input:
        val(batch)

    output:
        path(batch), emit: QCout_dir

    script:
    """
    Rscript ${moduleDir}/templates/QC.R --batch ${batch} \
        --INDIR "${params.INDIR}" --outdir "." --minSolidity ${params.minSolidity} \
        --minArea.um2 ${params.minArea.um2} --minnCount_RNA ${params.minnCount_RNA} \
        --doubletProportion ${params.doubletProportion}
    """
}

