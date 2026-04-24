process CREATE_PY_INPUTS {
    tag "${sample_id}_${rep_id}"

    label 'process_single_dynamic_memory'

    memory {
        def file_size = qced_rds.size()
        def mem = file_size <= 100_000_000 ? 15 :
                  file_size <= 200_000_000 ? 30 :
                  file_size <= 400_000_000 ? 60 :
                  file_size <= 800_000_000 ? 120 : 240
        "${mem * task.attempt} GB"
    }
    
    publishDir "${params.outdir}/qc_py_inputs/${sample_id}_${rep_id}", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), val(rep_id), path(qced_rds), path(qced_tf)
    
    output:
    tuple val(sample_id), val(rep_id), path(dir_py_inputs), emit: ch_py_inputs
    
    script:
    """
    ${projectDir}/scripts/create_py_inputs.R -s ${sample_id}_${rep_id} -r ${qced_rds} -c ${qced_tf} -o ${dir_py_inputs}
    """
}