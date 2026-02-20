process QC {
    publishDir "${params.outdir}/",
        mode: 'copy'
    tag "QC"
    label "lowMemST"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == "apptainer" ?
        'oras://community.wave.seqera.io/library/r-seurat-data_r-dplyr_r-ggplot2_r-patchwork_r-seurat:482e7aa7efc6ab6f':
        'community.wave.seqera.io/library/r-seurat-data_r-dplyr_r-ggplot2_r-patchwork_r-seurat:cdc203cd018c3c8b' }"

    input:
        val(batch)

    output:
        path(batch), emit: QCout_dir

    script:
    """
    Rscript ${moduleDir}/templates/QC.R --batch ${batch} --outdir "." --minSolidity ${params.minSolidity} --minArea.um2 ${params.minArea.um2} --minnCount_RNA ${params.minnCount_RNA}
    """
}

