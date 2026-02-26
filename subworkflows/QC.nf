include { QC } from '../modules/QC'
include {findDoublets} from '../modules/Doublets'

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
