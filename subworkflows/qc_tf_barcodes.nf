workflow qc_tf_barcodes {
    take:
    ch_tf_files

    main:
    QC_TF_BARCODES(ch_tf_files)
    ch_qced_stats = QC_TF_BARCODES.out.ch_qced_stats
    ch_qced_tf = QC_TF_BARCODES.out.ch_qced_tf

    FILTER_TF_BARCODES(ch_qced_tf)
    ch_filtered_tf = FILTER_TF_BARCODES.out.ch_filtered_tf
    if (params.top_n != 0) {
        ch_filtered_tf_top = FILTER_TF_BARCODES.out.ch_filtered_tf_top
    } else {
        ch_filtered_tf_top = Channel.empty()
    }

    emit:
    ch_qced_stats
    ch_qced_tf
    ch_filtered_tf
    ch_filtered_tf_top
}

process QC_TF_BARCODES {
    label 'process_high_dynamic_memory'

    memory {
        def file_size = read_1.size()
        def mem = file_size <= 10_000_000_000 ? 10 :
                  file_size <= 20_000_000_000 ? 20 :
                  file_size <= 40_000_000_000 ? 40 :
                  file_size <= 80_000_000_000 ? 80 : 160
        "${mem * task.attempt} GB"
    }

    publishDir "${params.outdir}/qc_tf/${sample_id}", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), val(run_id), path(read_1), path(read_2), path(barcode_csv), path(qced_cells)

    output:
    tuple val(sample_id), val(run_id), path("${sample_id}_${run_id}.tf_stats.tsv"), emit: ch_qced_stats
    tuple val(sample_id), val(run_id), path("${sample_id}_${run_id}.tf_barcodes.tsv"), emit: ch_qced_tf

    script:
    """
    python ${projectDir}/scripts/qc_tf_barcodes.py --reads          ${read_1},${read_2} \
                                                   --tf_barcode     ${barcode_csv} \
                                                   --cell_barcode   ${qced_cells} \
                                                   --tf_barcode_len ${params.tf_barcode_len} \
                                                   --marker_seq     ${params.marker_seq} \
                                                   --marker_start   ${params.marker_start} \
                                                   --marker_end     ${params.marker_end} \
                                                   --max_mismatch   ${params.max_mismatch} \
                                                   --output_prefix  ${sample_id}_${run_id}
    """
}

process FILTER_TF_BARCODES {
    label 'process_single_dynamic_memory'

    memory {
        def file_size = read_1.size()
        def mem = file_size <= 100_000_000 ? 4 :
                  file_size <= 1_000_000_000 ? 8 :
                  file_size <= 2_000_000_000 ? 16 :
                  file_size <= 4_000_000_000 ? 32 : 64
        "${mem * task.attempt} GB"
    }

    publishDir "${params.outdir}/qc_tf/${sample_id}", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), val(run_id), path(tf_barcode)

    output:
    tuple val(sample_id), val(run_id), path("${sample_id}_${run_id}.filtered_tfs.tsv"), emit: ch_filtered_tf
    tuple val(sample_id), val(run_id), path("${sample_id}_${run_id}.filtered_tfs_top.tsv"), emit: ch_filtered_tf_top, optional: true

    script:
    """
    ${projectDir}/scripts/filter_tf_barcodes.R -s ${sample_id}_${run_id} \
                                               -t ${tf_barcode} \
                                               --top_n ${params.top_n}
    """
}