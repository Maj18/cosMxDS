include { QC } from '../modules/QC'

workflow QC_workflow {
    take:
    batch

    main:
    QC(batch)
    ch_QCout_afterFiltering_dir = QC.out.QCout_afterFiltering_dir
    ch_QCout_beforeFiltering_dir = QC.out.QCout_beforeFiltering_dir
    ch_Seuratobject_filtered = QC.out.Seuratobject_filtered
    ch_Tissue_figs = QC.out.Tissue_figs
    ch_sessionInfo = QC.out.sessionInfo
    ch_QClog = QC.out.QClog
    
    emit:
    ch_QCout_afterFiltering_dir = ch_QCout_afterFiltering_dir
    ch_QCout_beforeFiltering_dir = ch_QCout_beforeFiltering_dir
    ch_Seuratobject_filtered = ch_Seuratobject_filtered
    ch_Tissue_figs = ch_Tissue_figs
    ch_sessionInfo = ch_sessionInfo
    ch_QClog = ch_QClog

}
