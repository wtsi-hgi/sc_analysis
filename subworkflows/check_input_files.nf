import java.util.zip.GZIPInputStream
import java.io.InputStreamReader
import java.io.BufferedReader

workflow check_input_files {
    take:
    ch_input

    main:
    CHECK_FILES(ch_input)
    ch_cr_files = CHECK_FILES.out.ch_cr_files
    ch_tf_files = CHECK_FILES.out.ch_tf_files

    emit:
    ch_cr_files
    ch_tf_files
}

process CHECK_FILES {
    label 'process_single'

    input:
    tuple val(sample_id), val(run_id), val(dir_cellrange_arc), val(r1_tf_barcodes), val(r2_tf_barcodes), val(tf_barcodes)

    output:
    tuple val(sample_id), path("${sample_id}.cr_gex.h5"), path("${sample_id}.cr_atac.tsv.gz"), emit: ch_cr_files
    tuple val(sample_id), path("${sample_id}.tf_barcodes.r1.fastq.gz"), path("${sample_id}.tf_barcodes.r2.fastq.gz"), path("${sample_id}.tf_barcodes.csv"), emit: ch_tf_files

    script:
    def file_cr_gex = file("${dir_cellrange_arc}/filtered_feature_bc_matrix.h5")
    def file_cr_atac = file("${dir_cellrange_arc}/atac_fragments.tsv.gz")
    def file_tf_read1 = file("${r1_tf_barcodes}")
    def file_tf_read2 = file("${r2_tf_barcodes}")
    def file_tf_barcodes = file("${tf_barcodes}")

    def valid_read_ext = [".fq", ".fastq", ".fq.gz", ".fastq.gz"]
    def valid_ref_ext = [".fa", ".fasta"]

    if (!file_cr_gex.exists()) {
        log.error("Error: ${file_cr_gex} is not found in ${dir_cellrange_arc}.")
        exit 1
    }

    if (!file_cr_atac.exists()) {
        log.error("Error: ${file_cr_atac} is not found in ${dir_cellrange_arc}.")
        exit 1
    }

    if (file_tf_read1.exists()) {
        if (!valid_read_ext.any { file_tf_read1.toString().endsWith(it) }) {
            log.error("Error: File format for ${file_tf_read1} is incorrect. Expected one of: ${valid_read_ext.join(', ')}")
            exit 1
        }
    } else {
        log.error("Error: ${file_tf_read1} is not found.")
        exit 1
    }

    if (file_tf_read2.exists()) {
        if (!valid_read_ext.any { file_tf_read2.toString().endsWith(it) }) {
            log.error("Error: File format for ${file_tf_read2} is incorrect. Expected one of: ${valid_read_ext.join(', ')}")
            exit 1
        }
    } else {
        log.error("Error: ${file_tf_read2} is not found.")
        exit 1
    }

    if (file_tf_barcodes.exists()) {
        def firstLine
        if (file_tf_barcodes.toString().endsWith(".gz")) {
            firstLine = new BufferedReader( new InputStreamReader(new GZIPInputStream(file_tf_barcodes.newInputStream()))).readLine()
        } else {
            firstLine = file_tf_barcodes.withReader { it.readLine() }
        }

        if (firstLine.contains(",")) {
            def header = firstLine.split(",").collect { it.toLowerCase() }
            if (header[0] != "barcode" || header[1] != "tf_name") {
                log.error("Error: ${file_tf_barcodes} file format is incorrect. Expected header: barcode,tf_name")
                exit 1
            }
        } else {
            log.error("Error: expect ${file_tf_barcodes} is a csv file.")
            exit 1
        }
    } else {
        log.error("Error: ${file_tf_barcodes} is not found.")
        exit 1
    }

    """
    echo "Checking: ${sample_id}"

    ln -s ${file_cr_gex} ${sample_id}.cr_gex.h5
    ln -s ${file_cr_atac} ${sample_id}.cr_atac.tsv.gz
    ln -s ${file_tf_read1} ${sample_id}.tf_barcodes.r1.fastq.gz
    ln -s ${file_tf_read2} ${sample_id}.tf_barcodes.r2.fastq.gz
    ln -s ${file_tf_barcodes} ${sample_id}.tf_barcodes.csv
    """
}
