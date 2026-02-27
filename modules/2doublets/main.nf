process findDoublets {
    publishDir "${params.outdir}/",
        mode: 'copy'
    tag "doublet"
    label "highMemMT"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == "apptainer" ?
        'docker://yuanli202004/seurat5.4_doubletfinder:latest':
        'docker://yuanli202004/seurat5.4_doubletfinder:latest' }"

    input:
        tuple(val(batch), file(INFILE))

    output:
        tuple(val(batch), path("${batch}/Doublets")), emit: QCdoublet_dir
        tuple(val(batch), path("${batch}/${batch}_filtered.RDS")), emit: filteredSeuratObject_wDoubletsInfo

    script:
    """
    Rscript ${moduleDir}/templates/findDoublets.R --batch ${batch} \
        --INFILE "${INFILE}" --outdir "." \
        --numCore ${params.cpus} \
        --doubletProportion ${params.doubletProportion}
    """
}

