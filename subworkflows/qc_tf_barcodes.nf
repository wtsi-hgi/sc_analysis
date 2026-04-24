include { QC_TF_BARCODES } from "$projectDir/modules/local/qc_tf_barcodes/main"
include { FILTER_TF_BARCODES } from "$projectDir/modules/local/filter_tf_barcodes/main"

workflow qc_tf_barcodes {
    take:
    ch_tf_files

    main:
    QC_TF_BARCODES(ch_tf_files)
    ch_qced_stats = QC_TF_BARCODES.out.ch_qced_stats
    ch_qced_tf = QC_TF_BARCODES.out.ch_qced_tf

    FILTER_TF_BARCODES(ch_qced_tf)
    ch_tf_cutoff_plots = FILTER_TF_BARCODES.out.ch_tf_cutoff_plots
    ch_filtered_tf = FILTER_TF_BARCODES.out.ch_filtered_tf
    ch_tf_filter_scatter = FILTER_TF_BARCODES.out.ch_tf_filter_scatter
    ch_tf_filter_boxplot = FILTER_TF_BARCODES.out.ch_tf_filter_boxplot
    if (params.top_n != 0) {
        ch_filtered_tf_top = FILTER_TF_BARCODES.out.ch_filtered_tf_top
    } else {
        ch_filtered_tf_top = Channel.empty()
    }

    emit:
    ch_qced_stats
    ch_qced_tf
    ch_tf_cutoff_plots
    ch_filtered_tf
    ch_tf_filter_scatter
    ch_tf_filter_boxplot
    ch_filtered_tf_top
}
