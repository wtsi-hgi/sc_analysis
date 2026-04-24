process GET_ATAC_DOUBLETS {
    tag "${sample_id}_${rep_id}"
    
    label 'process_single_dynamic_memory'

    memory {
        def file_size = file_atac.size()
        def mem = file_size <= 1_000_000_000 ? 2 :
                  file_size <= 2_000_000_000 ? 4 :
                  file_size <= 4_000_000_000 ? 8 :
                  file_size <= 8_000_000_000 ? 16 : 32
        "${mem * task.attempt} GB"
    }

    input:
    tuple val(sample_id), val(rep_id), path(file_raw), path(file_gex), path(file_atac), path(file_atac_tbi)

    output:
    tuple val(sample_id), val(rep_id), path("${sample_id}_${rep_id}.atac_doublets.tsv.gz"), emit: ch_atac_doublets

    script:
    """
    zcat ${file_atac} | grep -P '^chr([1-9]|1[0-9]|2[0-2])\t' | gzip > fragments_autosomes.tsv.gz

    bedtools intersect -a fragments_autosomes.tsv.gz \
                       -b ${projectDir}/data/hg38_repeatmasker_ucsc.tsv.gz \
                       -v | bgzip > fragments_clean.tsv.gz

    tabix -p bed fragments_clean.tsv.gz

    ${projectDir}/scripts/get_atac_doublets.R -s ${sample_id}_${rep_id} -f fragments_clean.tsv.gz

    gzip ${sample_id}_${rep_id}.atac_doublets.tsv
    """
}