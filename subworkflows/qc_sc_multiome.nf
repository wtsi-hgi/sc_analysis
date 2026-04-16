include { GET_ATAC_DOUBLETS } from "$projectDir/modules/local/get_atac_doublets/main"
include { QC_SC_MULTIOME } from "$projectDir/modules/local/qc_sc_multiome/main"

workflow qc_sc_multiome {
    take:
    ch_cr_files

    main:
    GET_ATAC_DOUBLETS(ch_cr_files)
    ch_atac_doublets = GET_ATAC_DOUBLETS.out.ch_atac_doublets

    ch_input = ch_cr_files.join(ch_atac_doublets, by: [0,1])

    QC_SC_MULTIOME(ch_input)
    ch_qced_object = QC_SC_MULTIOME.out.ch_qced_object
    ch_qced_summary = QC_SC_MULTIOME.out.ch_qced_summary
    ch_qced_cells = QC_SC_MULTIOME.out.ch_qced_cells

    emit:
    ch_qced_object
    ch_qced_summary
    ch_qced_cells
}
