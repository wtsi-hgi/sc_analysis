workflow generate_summary_report {
    take:
    ch_input

    main:
    CREATE_HTML_REPORT(ch_input)
    ch_qc_summary = CREATE_HTML_REPORT.out.ch_qc_summary

    emit:
    ch_qc_summary
}

process CREATE_HTML_REPORT {
    label 'process_single'

    publishDir "${params.outdir}/qc_report", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), val(rep_id), 
          val(qc_sc_stats), val(qc_sc_rna_plot), val(qc_sc_atac_plot), 
          val(qc_tf_stats), val(qc_tf_cutoff_plot), val(qc_tf_scatter_plot), val(qc_tf_boxplot_plot)

    output:
    tuple val(sample_id), path("${sample_id}_${rep_id}.qc_summary.html"), emit: ch_qc_summary

    script:
    """
    ${projectDir}/scripts/create_html_report.R -r ${projectDir}/scripts \
                                               -s ${sample_id}_${rep_id} \
                                               -t ${qc_sc_stats} \
                                               -m ${qc_sc_rna_plot} \
                                               -f ${qc_sc_atac_plot} \
                                               -a ${qc_tf_stats} \
                                               -c ${qc_tf_cutoff_plot} \
                                               -g ${qc_tf_scatter_plot} \
                                               -b ${qc_tf_boxplot_plot} \
                                               -w ${params.pipeline_name} \
                                               -v ${params.pipeline_version}

    """
}