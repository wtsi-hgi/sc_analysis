process QC_SC_MULTIOME {
    tag "${sample_id}_${rep_id}"
    
    label 'process_single_dynamic_memory'

    memory {
        def file_size = file_gex.size()
        def mem = file_size <= 100_000_000 ? 15 :
                  file_size <= 200_000_000 ? 30 :
                  file_size <= 400_000_000 ? 60 :
                  file_size <= 800_000_000 ? 120 : 240
        "${mem * task.attempt} GB"
    }

    publishDir "${params.outdir}/qc_sc/${sample_id}_${rep_id}", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), val(rep_id), path(file_raw), path(file_gex), path(file_atac), path(file_atac_tbi), path(file_atac_doublets)

    output:
    tuple val(sample_id), val(rep_id), path("${sample_id}_${rep_id}.qc_obj.rds"),  emit: ch_qced_object
    tuple val(sample_id), val(rep_id), path("${sample_id}_${rep_id}.qc_summary.tsv"), path("${sample_id}_${rep_id}.qc_rna_violin.png"), path("${sample_id}_${rep_id}.qc_atac_violin.png"), emit: ch_qced_summary
    tuple val(sample_id), val(rep_id), path("${sample_id}_${rep_id}.qc_cell_barcodes.tsv"), emit: ch_qced_cells

    script:
    def op_ambient = params.del_ambient ? "--del_ambient" : ""
    def op_doublet = params.mark_doublet ? "--mark_doublet --doublet_cells ${file_atac_doublets}" : ""
    """
    ${projectDir}/scripts/qc_sc_multiome.R -s ${sample_id}_${rep_id} \
                                           -r ${file_raw} \
                                           -g ${file_gex} \
                                           -a ${file_atac} ${op_ambient} ${op_doublet}
    """
}