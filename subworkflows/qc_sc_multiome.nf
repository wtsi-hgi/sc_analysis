workflow qc_sc_multiome {
    take:
    ch_cr_files

    main:
    QC_SC_MULTIOME(ch_cr_files)
    ch_qced_object = QC_SC_MULTIOME.out.ch_qced_object
    ch_qced_summary = QC_SC_MULTIOME.out.ch_qced_summary
    ch_qced_cells = QC_SC_MULTIOME.out.ch_qced_cells

    emit:
    ch_qced_object
    ch_qced_summary
    ch_qced_cells
}

process QC_SC_MULTIOME {
    label 'process_single_dynamic_memory'

    memory {
        def file_size = file_gex.size()
        def mem = file_size <= 100_000_000 ? 20 :
                  file_size <= 200_000_000 ? 40 :
                  file_size <= 400_000_000 ? 80 :
                  file_size <= 800_000_000 ? 160 : 320
        "${mem * task.attempt} GB"
    }

    publishDir "${params.outdir}/qc_sc/${sample_id}", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), val(run_id), path(file_gex), path(file_atac)

    output:
    tuple val(sample_id), val(run_id), path("${sample_id}_${run_id}.qc_obj.rds"),  emit: ch_qced_object
    tuple val(sample_id), val(run_id), path("${sample_id}_${run_id}.qc_summary.tsv"), path("${sample_id}_${run_id}.qc_rna_violin.png"), path("${sample_id}_${run_id}.qc_atac_violin.png"), emit: ch_qced_summary
    tuple val(sample_id), val(run_id), path("${sample_id}_${run_id}.qc_cell_barcodes.tsv"), emit: ch_qced_cells

    script:
    """
    ${projectDir}/scripts/qc_sc_per_sample.R -s ${sample_id}_${run_id} \
                                             -g ${file_gex} \
                                             -a ${file_atac}
    """
}