process QC_TF_BARCODES {
    tag "${sample_id}_${rep_id}"
    
    label 'process_high_dynamic_memory'

    memory {
        def file_size = read_1.size()
        def mem = file_size <= 10_000_000_000 ? 4 :
                  file_size <= 20_000_000_000 ? 8 :
                  file_size <= 40_000_000_000 ? 16 :
                  file_size <= 80_000_000_000 ? 32 : 64
        "${mem * task.attempt} GB"
    }

    publishDir "${params.outdir}/qc_tf/${sample_id}_${rep_id}", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), val(rep_id), path(read_1), path(read_2), path(barcode_csv), path(qced_cells)

    output:
    tuple val(sample_id), val(rep_id), path("${sample_id}_${rep_id}.tf_stats.tsv"), emit: ch_qced_stats
    tuple val(sample_id), val(rep_id), path("${sample_id}_${rep_id}.tf_barcodes.tsv"), emit: ch_qced_tf

    script:
    def op_has_umi = params.has_umi ? "--has_umi" : ""
    """
    python ${projectDir}/scripts/qc_tf_barcodes.py --reads          ${read_1},${read_2} \
                                                   --tf_barcode     ${barcode_csv} \
                                                   --cell_barcode   ${qced_cells} \
                                                   --tf_barcode_len ${params.tf_barcode_len} \
                                                   --marker_seq     ${params.marker_seq} \
                                                   --marker_start   ${params.marker_start} \
                                                   --marker_end     ${params.marker_end} \
                                                   --max_mismatch   ${params.max_mismatch} \
                                                   --output_prefix  ${sample_id}_${rep_id} ${op_has_umi}
    """
}