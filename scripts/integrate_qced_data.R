#!/usr/bin/env Rscript
quiet_library <- function(pkg) { suppressMessages(suppressWarnings(library(pkg, character.only = TRUE))) }
package_1 <- c("optparse", "tidyverse", "data.table", "ggplot2", "ggrepel")
package_2 <- c("Seurat", "Signac", "SeuratWrappers", "clusterProfiler", "GenomicRanges", "EnsDb.Hsapiens.v86")
packages <- c(package_1, package_2)
invisible(lapply(packages, quiet_library))

option_list <- list(make_option(c("-s", "--sample_ids"),  type = "character", help = "list of sample IDs",              default = NULL),
                    make_option(c("-r", "--rds_files"),   type = "character", help = "list of seurat object rds files", default = NULL),
                    make_option(c("-t", "--tf_files"),    type = "character", help = "list of tf barcode count files",  default = NULL),
                    make_option(c("-o", "--output_dir"),  type = "character", help = "output directory",                default = getwd()),
                    make_option(c("-p", "--prefix"),      type = "character", help = "output prefix",                   default = "integrated_qced"))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if(length(commandArgs(trailingOnly = TRUE)) == 0)
{
    print_help(opt_parser)
    quit(status = 1)
}

# -- check options -- #
if(is.null(opt$sample_ids)) stop("-s, list of sample IDs is required!", call. = FALSE)
if(is.null(opt$rds_files))  stop("-r, list of seurat object rds file is required!", call. = FALSE)
if(is.null(opt$tf_files))   stop("-t, list of tf barcode count file is required!", call. = FALSE)

# -- outputs -- #
if(!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)
setwd(opt$output_dir)

sample_prefix <- opt$prefix

#-- processing --#
sample_ids <- strsplit(opt$sample_ids, ",")[[1]]

message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "reading RDS ...")
list_obj_files <- strsplit(opt$rds_files, ",")[[1]]
list_objs <- lapply(list_obj_files, readRDS)

message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "reading tsv ...")
list_tf_files <- strsplit(opt$tf_files, ",")[[1]]
list_tf_counts <- lapply(list_tf_files, fread)

message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "integrating objects ...")
features <- SelectIntegrationFeatures(object.list = list_objs, nfeatures = 3000)
list_objs <- PrepSCTIntegration(object.list = list_objs, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = list_objs, normalization.method = "SCT", anchor.features = features)
integrated_obj <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

rm(list_objs)
invisible(gc(verbose = FALSE))

message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "integrating TF counts ...")
tf_counts <- rbindlist(list_tf_counts)
tf_mat <- dcast(tf_counts, tf_name ~ cell_barcode, value.var = "count", fill = 0)
tf_mat <- tf_mat[!is.na(tf_name)]
tf_mat <- as.data.frame(tf_mat)
rownames(tf_mat) <- tf_mat$tf_name
tf_mat$tf_name <- NULL
tf_mat <- as.matrix(tf_mat[, colnames(integrated_obj)])

rm(list_tf_counts)
rm(tf_counts)
invisible(gc(verbose = FALSE))

message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "creating final object ...")
tf_assay <- CreateAssayObject(data = tf_mat)
integrated_obj[["TF"]] <- tf_assay
integrated_obj[["TF"]]@counts <- tf_mat

rm(tf_mat)
rm(tf_assay)
invisible(gc(verbose = FALSE))

message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "creating output files ...")
saveRDS(integrated_obj, file = paste0(sample_prefix, ".integrated_qced.rds"))
