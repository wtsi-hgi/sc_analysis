#!/usr/bin/env Rscript
quiet_library <- function(pkg) { suppressMessages(suppressWarnings(library(pkg, character.only = TRUE))) }
packages <- c("optparse", "tidyverse", "data.table", "Seurat", "Signac", "SeuratWrappers", "SeuratDisk")
invisible(lapply(packages, quiet_library))

option_list <- list(make_option(c("-s", "--sample_id"),  type = "character", help = "sample ID",                   default = NULL),
                    make_option(c("-r", "--qced_rds"),   type = "character", help = "qced seurat object rds file", default = NULL),
                    make_option(c("-c", "--tf_count"),   type = "character", help = "tf barcode count file",       default = NULL),
                    make_option(c("-o", "--output_dir"), type = "character", help = "output directory",            default = getwd()),
                    make_option(c("-p", "--prefix"),     type = "character", help = "output prefix",               default = NULL))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if(length(commandArgs(trailingOnly = TRUE)) == 0)
{
    print_help(opt_parser)
    quit(status = 1)
}

# -- check options -- #
if(is.null(opt$sample_id)) stop("-s, sample ID is required!", call. = FALSE)
if(is.null(opt$qced_rds))  stop("-r, qced seurat object rds file is required!", call. = FALSE)
if(is.null(opt$tf_count))  stop("-c, tf barcode count file is required!", call. = FALSE)

# -- outputs -- #
if(!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)
setwd(opt$output_dir)

sample_prefix <- if (is.null(opt$prefix)) opt$sample_id else opt$prefix

# -- functions -- #
export_assay <- function(obj, assay, layer, outdir)
{
    if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

    DefaultAssay(obj) <- assay
    mat <- GetAssayData(obj, assay = assay, layer = layer)
    mat <- as(mat, "dgCMatrix")
  
    writeMM(mat, file.path(outdir, "matrix.mtx"))

    feature_type <- switch(assay,
                           "RNA"  = "Gene Expression",
                           "ATAC" = "Peaks",
                           "TF"   = "Counts",
                           assay)

    features <- data.frame(feature_id   = rownames(mat),
                           feature_name = rownames(mat),
                           feature_type = rep(feature_type, nrow(mat)))
    fwrite(features, file.path(outdir, "features.tsv"), sep = "\t", col.names = FALSE)

    cells <- data.frame(cell_id = colnames(mat))
    fwrite(cells, file.path(outdir, "barcodes.tsv"), sep = "\t", col.names = FALSE)
}

#-- processing --#
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "reading rds ...")
obj <- readRDS(opt$qced_rds)

message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "reading tsv ...")
tf_counts <- fread(opt$tf_count)

message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "integrating TF counts ...")
tf_mat <- dcast(tf_counts, tf_name ~ cell_barcode, value.var = "count", fill = 0)
tf_mat <- as.data.frame(tf_mat)
rownames(tf_mat) <- tf_mat$tf_name
tf_mat$tf_name <- NULL

common_cells <- intersect(colnames(tf_mat), colnames(obj))
obj <- subset(obj, cells = common_cells)
tf_mat <- as.matrix(tf_mat[, colnames(obj)])

message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "creating final object ...")
obj[["TF"]] <- CreateAssayObject(counts = tf_mat)

rm(tf_counts)
rm(tf_mat)
invisible(gc(verbose = FALSE))

message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "creating output files ...")
tmp_dir <- paste0(sample_prefix, "_py_inputs")
dir.create(tmp_dir, recursive = TRUE)

message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "    |--> exporting RNA assay ...")
export_assay(obj, "RNA", "counts", paste0(tmp_dir, "/", "rna"))

message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "    |--> exporting ATAC assay ...")
export_assay(obj, "ATAC", "counts", paste0(tmp_dir, "/", "atac"))

message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "    |--> exporting TF assay ...")
export_assay(obj, "TF", "counts", paste0(tmp_dir, "/", "tf"))
