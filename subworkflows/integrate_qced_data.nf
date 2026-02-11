workflow integrate_qced_data {
    take:
    ch_qced_files

    main:
    ch_input = ch_qced_files.groupTuple()
    INTEGRATE_QCED_DATA(ch_input)
    ch_integrated_qced = INTEGRATE_QCED_DATA.out.ch_integrated_qced

    emit:
    ch_integrated_qced
}

process INTEGRATE_QCED_DATA {
    label 'process_single_dynamic_memory'

    memory {
        def file_size_1 = qced_rds[0].size()
        def file_size_2 = qced_tf[0].size()
        def file_size_total = file_size_1 + file_size_2
        def mem = file_size_total <= 100_000_000 ? 16 :
                  file_size_total <= 1_000_000_000 ? 32 :
                  file_size_total <= 2_000_000_000 ? 64 :
                  file_size_total <= 4_000_000_000 ? 128 : 256
        "${mem * task.attempt} GB"
    }

    publishDir "${params.outdir}/qc_integration/${sample_id}", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(qced_rds), path(qced_tf)

    output:
    tuple val(sample_id), path("${sample_id}.integrated_qced.rds"), emit: ch_integrated_qced

    script:
    def list_run_ids = run_id.join(',')
    def list_qced_rds = qced_rds.join(',')
    def list_qced_tf = qced_tf.join(',')

    """
    ${projectDir}/scripts/integrate_qced_data.R -s ${list_run_ids} \
                                                -r ${list_qced_rds} \
                                                -t ${list_qced_tf} \
                                                -p ${sample_id}
    """



}