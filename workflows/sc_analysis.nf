/* ---- single cell data analysis pipeline ---- */

/* -- load modules -- */

/* -- load subworkflows -- */
include { check_input_files }   from '../subworkflows/check_input_files.nf'
include { qc_sc_multiome }      from '../subworkflows/qc_sc_multiome.nf'
include { qc_tf_barcodes }      from '../subworkflows/qc_tf_barcodes.nf'
include { integrate_qced_data } from '../subworkflows/integrate_qced_data.nf'

/* -- define functions -- */
def helpMessage() {
    log.info """
Usage:
    nextflow run sc_analysis/main.nf --sample_sheet "/path/of/sample/sheet"

    Mandatory arguments:
        --sample_sheet                path of the sample sheet
        --outdir                      the directory path of output results, default: the current directory

    Optional arguments:
    SC data QC:
        --del_ambient                 whether to remove ambient RNA, default: false
        --mark_doublet                whether to mark doublets, default: false
    
    TF barcode QC:
        --tf_barcode_len              length of TF barcode, default: 24
        --marker_seq                  marker sequence for TF barcode, default: "GAAAGGACGA"
        --marker_start                start position of marker sequence, default: 25
        --marker_end                  end position of marker sequence, default: 50
        --max_mismatch                maximum mismatch allowed for match sequence, default: 1
        --top_n                       keep top N TFs if 0 keep all, default: 10

    """
}

/* -- initialising parameters -- */
params.help           = null

params.sample_sheet   = null
params.outdir         = params.outdir.        ?: "$PWD"

params.del_ambient    = params.del_ambient    ?: false
params.mark_doublet   = params.mark_doublet   ?: false

params.tf_barcode_len = params.tf_barcode_len ?: 24
params.marker_seq     = params.marker_seq     ?: "GAAAGGACGA"
params.marker_start   = params.marker_start   ?: 25
params.marker_end     = params.marker_end     ?: 50
params.max_mismatch   = params.max_mismatch   ?: 1
params.top_n          = params.top_n          ?: null

/* -- check parameters -- */
if (params.help) {
    helpMessage()
    exit 0
}

if (params.sample_sheet) {
    def sep = params.sample_sheet.endsWith('.tsv') ? '\t' : ','
    ch_input = Channel.fromPath(file(params.sample_sheet), checkIfExists: true)
                      .splitCsv(header: true, sep: sep)
    
    def required_cols = ['sample_id', 'run_id', 'dir_cellrange_arc', 'r1_tf_barcodes', 'r2_tf_barcodes', 'tf_barcodes']
    def header_line = new File(params.sample_sheet).readLines().head()
    def header = header_line.split(sep)
    def missing = required_cols.findAll { !(it in header) }

    if (missing) {
        log.error "Sample sheet is missing required columns: ${missing.join(', ')}"
        exit 1
    }
} else {
    helpMessage()
    log.info("Error: Please specify the full path of the sample sheet!\n")
    exit 1
}

if (!file(params.outdir).isDirectory()) {
    log.error("Invalid output directory: ${params.outdir}. Please specify a valid directory.")
    exit 1
}

/* -- workflow -- */
workflow sc_analysis {
    /* -- check inputs -- */
    check_input_files(ch_input)
    ch_cr_files = check_input_files.out.ch_cr_files
    ch_tf_files = check_input_files.out.ch_tf_files

    /* -- step 1: QC single cell multiome data -- */
    qc_sc_multiome(ch_cr_files)
    ch_qced_object = qc_sc_multiome.out.ch_qced_object
    ch_qced_cells = qc_sc_multiome.out.ch_qced_cells

    /* -- step 2: QC TF barcodes -- */
    ch_input = ch_tf_files.join(ch_qced_cells, by: [0,1])
    qc_tf_barcodes(ch_input)
    ch_qced_tf = params.top_n == 0 ? qc_tf_barcodes.out.ch_filtered_tf : qc_tf_barcodes.out.ch_filtered_tf_top

    /* -- step 3: integration -- */
    ch_input = ch_qced_object.join(ch_qced_tf, by: [0,1])
    integrate_qced_data(ch_input)
    ch_integrated_qced = integrate_qced_data.out.ch_integrated_qced
}