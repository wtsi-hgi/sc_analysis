workflow qc_sc_multiome {
    take:
    ch_cr_files

    main:
    GET_ATAC_DOUBLETS(ch_cr_files)
    ch_atac_doublets = GET_ATAC_DOUBLETS.out.ch_atac_doublets

    ch_input = ch_cr_files.join(ch_atac_doublets, by: [0,1])

    QC_SC_MULTIOME(ch_input)
    ch_qced_object = QC_SC_MULTIOME.out.ch_qced_object
    ch_qced_summary = QC_SC_MULTIOME.out.ch_qced_summary
    ch_qced_cells = QC_SC_MULTIOME.out.ch_qced_cells

    emit:
    ch_qced_object
    ch_qced_summary
    ch_qced_cells
}

process GET_ATAC_DOUBLETS {
    label 'process_single_dynamic_memory'

    memory {
        def file_size = file_atac.size()
        def mem = file_size <= 100_000_000 ? 20 :
                  file_size <= 200_000_000 ? 40 :
                  file_size <= 400_000_000 ? 80 :
                  file_size <= 800_000_000 ? 160 : 320
        "${mem * task.attempt} GB"
    }

    input:
    tuple val(sample_id), val(run_id), path(file_raw), path(file_gex), path(file_atac), path(file_atac_tbi)

    output:
    tuple val(sample_id), val(run_id), path("${sample_id}_${run_id}.atac_doublets.tsv.gz"), emit: ch_atac_doublets

    script:
    """
    zcat ${file_atac} | grep -P '^chr([1-9]|1[0-9]|2[0-2])\t' | gzip > fragments_autosomes.tsv.gz

    bedtools intersect -a fragments_autosomes.tsv.gz \
                       -b ${projectDir}/data/hg38_repeatmasker_ucsc.tsv.gz \
                       -v | bgzip > fragments_clean.tsv.gz

    tabix -p bed fragments_clean.tsv.gz

    ${projectDir}/scripts/get_atac_doublets.R -s ${sample_id}_${run_id} -f fragments_clean.tsv.gz

    gzip ${sample_id}_${run_id}.atac_doublets.tsv
    """
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
    tuple val(sample_id), val(run_id), path(file_raw), path(file_gex), path(file_atac), path(file_atac_tbi), path(file_atac_doublets)

    output:
    tuple val(sample_id), val(run_id), path("${sample_id}_${run_id}.qc_obj.rds"),  emit: ch_qced_object
    tuple val(sample_id), val(run_id), path("${sample_id}_${run_id}.qc_summary.tsv"), path("${sample_id}_${run_id}.qc_rna_violin.png"), path("${sample_id}_${run_id}.qc_atac_violin.png"), emit: ch_qced_summary
    tuple val(sample_id), val(run_id), path("${sample_id}_${run_id}.qc_cell_barcodes.tsv"), emit: ch_qced_cells

    script:
    def op_ambient = params.del_ambient ? "--del_ambient" : ""
    def op_doublet = params.mark_doublet ? "--mark_doublet --doublet_cells ${file_atac_doublets}" : ""
    """
    ${projectDir}/scripts/qc_sc_per_sample.R -s ${sample_id}_${run_id} \
                                             -r ${file_raw} \
                                             -g ${file_gex} \
                                             -a ${file_atac} ${op_ambient} ${op_doublet}
    """
}