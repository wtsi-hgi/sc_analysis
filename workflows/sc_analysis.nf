/* ---- single cell data analysis pipeline ---- */

/* -- load modules -- */

/* -- load subworkflows -- */
include { check_input_files }         from '../subworkflows/check_input_files.nf'

/* -- define functions -- */
def helpMessage() {
    log.info """
Usage:
    nextflow run sc_analysis/main.nf --sample_sheet "/path/of/sample/sheet"

    Mandatory arguments:
        --sample_sheet                path of the sample sheet
        --outdir                      the directory path of output results, default: the current directory
    """
}

/* -- initialising parameters -- */
params.help                        = null

params.sample_sheet                = null
params.outdir                      = params.outdir                      ?: "$PWD"

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

    /* -- step 1: QC TF barcodes -- */


    /* -- step 2: QC single cell multiome data -- */


    /* -- step 3: integration -- */

}