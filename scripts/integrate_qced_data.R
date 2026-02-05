#!/usr/bin/env Rscript
quiet_library <- function(pkg) { suppressMessages(suppressWarnings(library(pkg, character.only = TRUE))) }
package_1 <- c("optparse", "tidyverse", "data.table", "ggplot2", "ggrepel")
package_2 <- c("Seurat", "Signac", "SeuratWrappers", "clusterProfiler", "GenomicRanges", "EnsDb.Hsapiens.v86")
packages <- c(package_1, package_2)
invisible(lapply(packages, quiet_library))

option_list <- list(make_option(c("-s", "--sample_id"),  type = "character", help = "sample ID",                       default = NULL),
                    make_option(c("-r", "--rds_files"),  type = "character", help = "list of seurat object rds files", default = NULL),
                    make_option(c("-t", "--tf_files"),   type = "character", help = "list of tf barcode count files",  default = NULL),
                    make_option(c("-o", "--output_dir"), type = "character", help = "output directory",                default = getwd()),
                    make_option(c("-p", "--prefix"),     type = "character", help = "output prefix",                   default = NULL))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if(length(commandArgs(trailingOnly = TRUE)) == 0)
{
    print_help(opt_parser)
    quit(status = 1)
}

# -- check options -- #
if(is.null(opt$sample_id)) stop("-s, sample ID is required!", call. = FALSE)
if(is.null(opt$rds_file))  stop("-r, list of seurat object rds file is required!", call. = FALSE)
if(is.null(opt$tf_file))   stop("-t, list of tf barcode count file is required!", call. = FALSE)

# -- outputs -- #
if(!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)
setwd(opt$output_dir)

sample_prefix <- ifelse(is.null(opt$prefix), opt$sample_id, opt$prefix)

#-- processing --#
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "Sample QC: ", opt$sample_id, " --> reading inputs ...")
