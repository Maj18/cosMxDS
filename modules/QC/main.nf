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
        path("${batch}/afterFiltering"), emit: QCout_afterFiltering_dir
        path("${batch}/beforeFiltering"), emit: QCout_beforeFiltering_dir
        path("${batch}/Eleni_Male_filtered.RDS"), emit: Seuratobject_filtered
        path("${batch}/Tissue.pdf"), emit: Tissue_figs
        path("${batch}/sessionInfo.txt"), emit: sessionInfo
        path("${batch}/QClog.txt"), emit: QClog

    script:
    """
    Rscript ${moduleDir}/templates/QC.R --batch ${batch} \
        --INDIR "${params.INDIR}" --outdir "." --minSolidity ${params.minSolidity} \
        --minArea.um2 ${params.minArea.um2} --minnCount_RNA ${params.minnCount_RNA} \
        --numCore ${params.cpus} \
        --doubletProportion ${params.doubletProportion}
    """
}

