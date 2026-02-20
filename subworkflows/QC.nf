include { QC } from '../modules/QC'

workflow QC_workflow {
    take:
    batch

    main:
    QC(batch)
    ch_QCout = QC.out.QCout_dir
    
    emit:
    ch_QCout = ch_QCout
}
