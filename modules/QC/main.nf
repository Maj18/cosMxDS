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
        tuple(val(batch), path("${batch}/afterFiltering")), emit: QCout_afterFiltering_dir
        tuple(val(batch), path("${batch}/beforeFiltering")), emit: QCout_beforeFiltering_dir
        tuple(val(batch), path("${batch}/Eleni_Male_filtered.RDS")), emit: SeuratObject_filtered
        tuple(val(batch), path("${batch}/Tissue.pdf")), emit: Tissue_figs
        tuple(val(batch), path("${batch}/sessionInfo.txt")), emit: sessionInfo
        tuple(val(batch), path("${batch}/QClog.txt")), emit: QClog

    script:
    """
    Rscript ${moduleDir}/templates/QC.R --batch ${batch} \
        --INDIR "${params.INDIR}" --outdir "." --minSolidity ${params.minSolidity} \
        --minArea.um2 ${params.minArea.um2} --minnCount_RNA ${params.minnCount_RNA}
    """
}

