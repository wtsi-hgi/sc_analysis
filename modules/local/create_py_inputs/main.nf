process CREATE_PY_INPUTS {
    tag "${sample_id}_${rep_id}"

    label 'process_single_dynamic_memory'

    memory {
        def file_size = qced_rds.size()
        def mem = file_size <= 1_000_000_000 ? 20 :
                  file_size <= 2_000_000_000 ? 40 :
                  file_size <= 4_000_000_000 ? 80 :
                  file_size <= 8_000_000_000 ? 160 : 320
        "${mem * task.attempt} GB"
    }
    
    publishDir "${params.outdir}/qc_py_inputs/${sample_id}_${rep_id}", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), val(rep_id), path(qced_rds), path(qced_tf)
    
    output:
    tuple val(sample_id), val(rep_id), path("${sample_id}_${rep_id}_py_inputs"), emit: ch_py_inputs
    
    script:
    """
    ${projectDir}/scripts/create_py_inputs.R -s ${sample_id}_${rep_id} -r ${qced_rds} -c ${qced_tf}

    gzip ${sample_id}_${rep_id}_py_inputs/*/*.tsv
    gzip ${sample_id}_${rep_id}_py_inputs/*/*.mtx
    """
}