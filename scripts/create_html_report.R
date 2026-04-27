#!/usr/bin/env Rscript
quiet_library <- function(pkg) { suppressMessages(suppressWarnings(library(pkg, character.only = TRUE))) }
packages <- c("tidyverse", "data.table", "vroom", "ggVennDiagram", "htmltools", "reactable", "optparse", "sparkline", "UpSetR", "patchwork", "glue", "scales", "ggExtra", "gtools")
invisible(lapply(packages, quiet_library))

# -- options -- #
option_list <- list(make_option(c("-r", "--rscript_dir"),        type = "character", help = "directory path of R scripts",     default = NULL),
                    make_option(c("-s", "--sample_id"),          type = "character", help = "sample ID",                       default = NULL),
                    make_option(c("-t", "--qc_sc_stats"),        type = "character", help = "QC single-cell stats file",       default = NULL),
                    make_option(c("-m", "--qc_sc_rna_plot"),     type = "character", help = "QC single-cell RNA plots files",  default = NULL),
                    make_option(c("-f", "--qc_sc_atac_plot"),    type = "character", help = "QC single-cell ATAC plots files", default = NULL),
                    make_option(c("-a", "--qc_tf_stats"),        type = "character", help = "QC TF stats files",               default = NULL),
                    make_option(c("-c", "--qc_tf_cutoff_plot"),  type = "character", help = "QC TF cutoff plots files",        default = NULL),
                    make_option(c("-g", "--qc_tf_scatter_plot"), type = "character", help = "QC TF scatter plots files",      default = NULL),
                    make_option(c("-b", "--qc_tf_boxplot_plot"), type = "character", help = "QC TF boxplot plots files",     default = NULL),
                    make_option(c("-o", "--output_dir"),         type = "character", help = "output directory",                default = getwd()),
                    make_option(c("-p", "--prefix"),             type = "character", help = "output prefix",                   default = NULL),
                    make_option(c("-w", "--pl_name"),            type = "character", help = "pipeline name",                   default = "sc_analysis"),
                    make_option(c("-v", "--pl_version"),         type = "character", help = "pipeline version",                default = "dev"))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if(length(commandArgs(trailingOnly = TRUE)) == 0)
{
    print_help(opt_parser)
    quit(status = 1)
}

# -- check options -- #
if(is.null(opt$rscript_dir))        stop("-r, directory path of R scripts is required!", call. = FALSE)
if(is.null(opt$sample_id))          stop("-s, sample ID is required!", call. = FALSE)
if(is.null(opt$qc_sc_stats))        stop("-t, QC single-cell stats file is required!", call. = FALSE)
if(is.null(opt$qc_sc_rna_plot))    stop("-m, QC single-cell RNA plots files is required!", call. = FALSE)
if(is.null(opt$qc_sc_atac_plot))   stop("-f, QC single-cell ATAC plots files is required!", call. = FALSE)
if(is.null(opt$qc_tf_stats))        stop("-a, QC TF stats files is required!", call. = FALSE)
if(is.null(opt$qc_tf_cutoff_plot)) stop("-c, QC TF cutoff plots files is required!", call. = FALSE)
if(is.null(opt$qc_tf_scatter_plot)) stop("-g, QC TF scatter plots files is required!", call. = FALSE)
if(is.null(opt$qc_tf_boxplot_plot)) stop("-b, QC TF boxplot plots files is required!", call. = FALSE)

# -- modules -- #
source(file.path(opt$rscript_dir, "report_utils.R"))
source(file.path(opt$rscript_dir, "report_html.R"))

# -- inputs -- #

# -- outputs -- #
if(!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)
setwd(opt$output_dir)

sample_prefix <- if(!is.null(opt$prefix)) opt$prefix else opt$sample_id

# -- reporting -- #
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "Creating final html report...")

file_render_context <- paste0(sample_prefix, ".qced_summary.Rmd")
create_html_render(opt$pl_name,
                   opt$pl_version,
                   opt$qc_sc_stats, 
                   opt$qc_sc_rna_plot,
                   opt$qc_sc_atac_plot,
                   opt$qc_tf_stats,
                   opt$qc_tf_cutoff_plot,
                   opt$qc_tf_scatter_plot,
                   opt$qc_tf_boxplot_plot,
                   file_render_context)

rmarkdown::render(file_render_context, clean = TRUE, quiet = TRUE)
invisible(file.remove(file_render_context))
