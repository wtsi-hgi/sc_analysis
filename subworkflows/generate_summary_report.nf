workflow generate_summary_report {
    take:
    ch_input

    main:
    ch_sample.map { sample_id, run_id, qc_sc_stats, qc_sc_rna_plots, qc_sc_atac_plots, qc_tf_stats, qc_tf_cutoff_plots }
}

process CREATE_HTML_REPORT {
    label 'process_single_dynamic_memory'

    memory {
        def file_size = qc_sc_stats.size() + qc_tf_stats.size()
        def mem = file_size <= 10_000_000_000 ? 4 :
                  file_size <= 20_000_000_000 ? 8 :
                  file_size <= 40_000_000_000 ? 16 :
                  file_size <= 80_000_000_000 ? 32 : 64
        "${mem * task.attempt} GB"
    }

    publishDir "${params.outdir}/qc_report", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), val(run_id), val(qc_sc_stats), val(qc_sc_rna_plots), val(qc_sc_atac_plots), val(qc_tf_stats), val(qc_tf_cutoff_plots)

    output:
    tuple val(sample_id), path("${sample_id}.qced_summary.html"), emit: ch_qced_summary

    script:
    def list_run_ids = run_id.join(',')
    def list_qc_sc_stats = qc_sc_stats.join(',')
    def list_qc_sc_rna_plots = qc_sc_rna_plots.join(',')
    def list_qc_sc_atac_plots = qc_sc_atac_plots.join(',')
    def list_qc_tf_stats = qc_tf_stats.join(',')
    def list_qc_tf_cutoff_plots = qc_tf_cutoff_plots.join(',')

    """

    """
}