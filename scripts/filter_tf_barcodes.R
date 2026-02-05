#!/usr/bin/env Rscript
quiet_library <- function(pkg) { suppressMessages(suppressWarnings(library(pkg, character.only = TRUE))) }
package_1 <- c("optparse", "tidyverse", "data.table", "ggplot2")
package_2 <- c("reshape2", "scales", "ggExtra", "ggrepel")
packages <- c(package_1, package_2)
invisible(lapply(packages, quiet_library))

option_list <- list(make_option(c("-s", "--sample_id"),  type = "character", help = "sample ID",                     default = NULL),
                    make_option(c("-t", "--tf_barcode"), type = "character", help = "TF barcode file",               default = NULL),
                    make_option(c("-o", "--output_dir"), type = "character", help = "output directory",              default = getwd()),
                    make_option(c("-p", "--prefix"),     type = "character", help = "output prefix",                 default = NULL),
                    make_option("--top_n",               type = "integer",   help = "keep top N TFs, if 0 keep all", default = 0))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if(length(commandArgs(trailingOnly = TRUE)) == 0)
{
    print_help(opt_parser)
    quit(status = 1)
}

# -- check options -- #
if(is.null(opt$sample_id))  stop("-s, sample ID is required!", call. = FALSE)
if(is.null(opt$tf_barcode)) stop("-t, TF barcode file is required!", call. = FALSE)

# -- inputs -- #
tf_counts <- fread(opt$tf_barcode)

# -- outputs -- #
if(!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)
setwd(opt$output_dir)

sample_prefix <- ifelse(is.null(opt$prefix), opt$sample_id, opt$prefix)

#-- processing --#
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "Sample QC: ", opt$sample_id, " --> keeping TF barcodes with counts > 1 ...")
tf_counts <- tf_counts[, .(count = sum(count)), by = .(cell_barcode, tf_name)]
tf_counts <- tf_counts[count > 1]

message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "Sample QC: ", opt$sample_id, " --> clustering TF barcodes by the counts ...")
tf_counts[, cluster := {
    if (.N > 10) {
        ck <- Ckmeans.1d.dp(log2(count), k = 2, y = 1)$cluster
        ck
    } else {
        rep(2L, .N)
    }
}, by = cell_barcode]

message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "Sample QC: ", opt$sample_id, " --> filtering TF barcodes by the clusters ...")
tf_counts_f <- tf_counts[cluster == 2]
tf_counts_f[, cluster := NULL]

fwrite(tf_counts_f, paste0(sample_prefix, ".filtered_tfs.tsv"), row.names = F, sep = "\t")

tf_counts_f_n_tf <- tf_counts_f[, .(n_tf         = uniqueN(tf_name), 
                                    mean_count   = as.numeric(mean(count)), 
                                    median_count = as.numeric(median(count))), by = .(cell_barcode)]

png(paste0(sample_prefix, ".filtered_tfs.scatter.png"), width = 2000, height = 1600, res = 150)
p <- ggplot(tf_counts_f_n_tf, aes(n_tf, mean_count)) +
        geom_point(shape = 19, alpha = 0.4, size = 0.8, col = "royalblue") +
        theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
        theme(axis.title = element_text(size = 20, face = "bold", family = "Arial")) +
        theme(plot.title = element_text(size = 20, face = "bold.italic", family = "Arial")) +
        theme(axis.text = element_text(size = 16, face = "bold")) +
        scale_x_continuous(breaks = replace(seq(0, ceiling(max(tf_counts_f_n_tf$n_tf)), by = 30), 1, 10)) +
        scale_y_log10() +
        geom_vline(xintercept = 10, color = "tomato", linetype = "dashed", linewidth = 0.8)
ggMarginal(p, type = "density", col = "brown1", fill = "brown1", alpha = 0.5)
dev.off()

bin_values <- c(seq(0, 50, 5), seq(60, 100 , 10), Inf)
bin_labels <- c("0~5",     "5~10",    "10~15",   "15~20", 
                "20~25",   "25~30",   "30~35",   "35~40",
                "40~45",   "45~50",   "50~60",   "60~70",
                "70~80",   "80~90",   "90~100",  "100+")
tf_counts_f_n_tf[, tf_bin := cut(n_tf, breaks = bin_values, labels = bin_labels, right=FALSE)]

tf_counts_f_boxplot <- melt(tf_counts_f_n_tf, 
                            id.vars = c("cell_barcode", "tf_bin"), 
                            measure.vars = c("mean_count", "median_count"),
                            variable.name = "count_type", value.name = "count_value")

png(paste0(sample_prefix, ".filtered_tfs.boxplot.png"), width = 2000, height = 1600, res = 150)
ggplot(tf_counts_f_boxplot, aes(x = tf_bin, y = count_value, fill = count_type)) +
    geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.8), outlier.colour = "darkgrey", linewidth = 0.1) +
    scale_y_log10() +
    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
    theme(axis.title = element_text(size = 10, face = "bold", family = "Arial")) +
    theme(plot.title = element_text(size = 10, face = "bold.italic", family = "Arial")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "TF Bin", y = "Count", title = "Count by TF Bin")
dev.off()

if(opt$top_n > 0)
{
    message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "Sample QC: ", opt$sample_id, " --> filtering TF barcodes by the clusters ...")
    tf_counts_top <- tf_counts_f[, {if (.N > opt$top_n) .SD[order(-count)][1 : opt$top_n] else .SD}, by = cell_barcode]

    fwrite(tf_counts_top, paste0(sample_prefix, ".filtered_tfs_top.tsv"), row.names = F, sep = "\t")

    tf_counts_top_n_tf <- tf_counts_top[, .(n_tf         = uniqueN(tf_name), 
                                            mean_count   = as.numeric(mean(count)), 
                                            median_count = as.numeric(median(count))), by = .(cell_barcode)]

    png(filename = paste0(sample_prefix, ".filtered_tfs_top.scatter.png"), width = 2000, height = 1600, res = 150)
    p <- ggplot(tf_counts_top_n_tf, aes(n_tf, mean_count)) +
            geom_point(shape = 19, alpha = 0.4, size = 0.8, col = "royalblue") +
            theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
            theme(axis.title = element_text(size = 20, face = "bold", family = "Arial")) +
            theme(plot.title = element_text(size = 20, face = "bold.italic", family = "Arial")) +
            theme(axis.text = element_text(size = 16, face = "bold")) +
            scale_x_continuous(breaks = replace(seq(0, ceiling(max(tf_counts_top_n_tf$n_tf)), by = 30), 1, 10)) +
            scale_y_log10() +
            geom_vline(xintercept = 10, color = "tomato", linetype = "dashed", linewidth = 0.8)
    ggMarginal(p, type = "density", col = "brown1", fill = "brown1", alpha = 0.5)
    dev.off()

    tf_counts_top_boxplot <- melt(tf_counts_top_n_tf, 
                                  id.vars = c("cell_barcode", "n_tf"), 
                                  measure.vars = c("mean_count", "median_count"),
                                  variable.name = "count_type", value.name = "count_value")

    png(paste0(sample_prefix, ".filtered_tfs_top.boxplot.png"), width = 2000, height = 1600, res = 150)                              
    ggplot(tf_counts_top_boxplot, aes(x = factor(n_tf), y = count_value, fill = count_type)) +
        geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.8), outlier.colour = "darkgrey", linewidth = 0.1) +
        scale_y_log10() +
        theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
        theme(axis.title = element_text(size = 10, face = "bold", family = "Arial")) +
        theme(plot.title = element_text(size = 10, face = "bold.italic", family = "Arial")) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(x = "TF Bin", y = "Count", title = "Count by TF Bin")
    dev.off()
}
