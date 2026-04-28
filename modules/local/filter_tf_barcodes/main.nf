process FILTER_TF_BARCODES {
    tag "${sample_id}_${rep_id}"

    label 'process_single_dynamic_memory'

    memory {
        def file_size = tf_barcode.size()
        def mem = file_size <= 100_000_000 ? 1 :
                  file_size <= 1_000_000_000 ? 2 :
                  file_size <= 2_000_000_000 ? 4 :
                  file_size <= 4_000_000_000 ? 8 : 16
        "${mem * task.attempt} GB"
    }

    publishDir "${params.outdir}/qc_tf/${sample_id}_${rep_id}", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), val(rep_id), path(tf_barcode)

    output:
    tuple val(sample_id), val(rep_id), path("${sample_id}_${rep_id}.tf_cutoff_plots.png"), emit: ch_tf_cutoff_plots
    tuple val(sample_id), val(rep_id), path("${sample_id}_${rep_id}.filtered_tfs.tsv"), emit: ch_filtered_tf
    tuple val(sample_id), val(rep_id), path("${sample_id}_${rep_id}.filtered_tfs_top.tsv"), emit: ch_filtered_tf_top, optional: true
    tuple val(sample_id), val(rep_id), path("${sample_id}_${rep_id}.filtered_tfs.scatter.png"), emit: ch_tf_filter_scatter
    tuple val(sample_id), val(rep_id), path("${sample_id}_${rep_id}.filtered_tfs.boxplot.png"), emit: ch_tf_filter_boxplot

    script:
    def op_has_umi = params.has_umi ? "--has_umi" : ""
    """
    ${projectDir}/scripts/filter_tf_barcodes.R -s ${sample_id}_${rep_id} \
                                               -t ${tf_barcode} \
                                               -c ${params.tf_cutoff} \
                                               --top_n ${params.top_n} ${op_has_umi}
    """
}