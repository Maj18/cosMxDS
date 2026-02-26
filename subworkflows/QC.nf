include { QC } from '../modules/QC'
include { findDoublets } from '../modules/Doublets'
include { mergeBatch } from '../modules/mergeBatch'

workflow QC_workflow {
    take:
    batch

    main:
    QC(batch)
    ch_QCout_afterFiltering_dir = QC.out.QCout_afterFiltering_dir
    ch_QCout_beforeFiltering_dir = QC.out.QCout_beforeFiltering_dir
    ch_SeuratObject_filtered = QC.out.SeuratObject_filtered
    ch_Tissue_figs = QC.out.Tissue_figs
    ch_sessionInfo = QC.out.sessionInfo
    ch_QClog = QC.out.QClog

    findDoublets(ch_SeuratObject_filtered)
    ch_QCdoublet_dir = findDoublets.out.QCdoublet_dir
    ch_filteredSeuratObject_wDoubletsInfo = findDoublets.out.filteredSeuratObject_wDoubletsInfo

    ch_combined = ch_filteredSeuratObject_wDoubletsInfo.collect().map { tuples ->
        def batch_str = tuples.collect { it[0] }.join(",")
        def file_str = tuples.collect { it[1] }.join(",")
        tuple(batch_str, file_str)
    }
    mergeBatch(ch_combined)

    
    emit:
    ch_QCout_afterFiltering_dir = ch_QCout_afterFiltering_dir
    ch_QCout_beforeFiltering_dir = ch_QCout_beforeFiltering_dir
    ch_Seuratobject_filtered = ch_Seuratobject_filtered
    ch_Tissue_figs = ch_Tissue_figs
    ch_sessionInfo = ch_sessionInfo
    ch_QClog = ch_QClog
    ch_QCdoublet_dir = ch_QCdoublet_dir
    ch_filteredSeuratObject_wDoubletsInfo = ch_filteredSeuratObject_wDoubletsInfo

}
