#!/usr/bin/env Rscript
quiet_library <- function(pkg) { suppressMessages(suppressWarnings(library(pkg, character.only = TRUE))) }
packages <- c("optparse", "tidyverse", "data.table", "scDblFinder")
invisible(lapply(packages, quiet_library))

option_list <- list(make_option(c("-s", "--sample_id"),  type = "character", help = "sample ID",               default = NULL),
                    make_option(c("-f", "--fragments"),  type = "character", help = "scATAC fragment gz file", default = NULL),
                    make_option(c("-o", "--output_dir"), type = "character", help = "output directory",        default = getwd()),
                    make_option(c("-p", "--prefix"),     type = "character", help = "output prefix",           default = NULL))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if(length(commandArgs(trailingOnly = TRUE)) == 0)
{
    print_help(opt_parser)
    quit(status = 1)
}

# -- check options -- #
if(is.null(opt$sample_id))  stop("-s, sample ID is required!", call. = FALSE)
if(is.null(opt$fragments))  stop("-f, fragments file is required!", call. = FALSE)

# -- outputs -- #
if(!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)
setwd(opt$output_dir)

sample_prefix <- ifelse(is.null(opt$prefix), opt$sample_id, opt$prefix)

#-- processing --#
res <- amulet(opt$fragments)
doublet_cells <- data.table(barcode = rownames(res[res$q.value < 0.01,]))

output_file <- paste0(sample_prefix, ".atac_doublets.tsv")
fwrite(doublet_cells, output_file, sep = "\t", col.names = FALSE)
