library(data.table)
library(ggplot2)
library(ggExtra)

input_dir <- "/lustre/scratch127/gengen/projects_v2/dual_tf/MORF_10/Multiome/analysis_by_hgi/2_tf_barcodes"

tf_barcodes <- fread(file.path(input_dir, "morf10_tf.bc_umi_per_cell.tsv.gz"))
n_tf_per_cell <- tf_barcodes[, .(n_tf       = uniqueN(tf_name), 
                                 mean_umi   = as.numeric(mean(umi_count)), 
                                 median_umi = as.numeric(median(umi_count))), by = .(cell_barcode)]

p <- ggplot(n_tf_per_cell, aes(n_tf, mean_umi)) +
        geom_point(shape = 19, alpha = 0.4, size = 0.8, col = "royalblue") +
        theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
        theme(axis.title = element_text(size = 20, face = "bold", family = "Arial")) +
        theme(plot.title = element_text(size = 20, face = "bold.italic", family = "Arial")) +
        theme(axis.text = element_text(size = 16, face = "bold")) +
        scale_x_continuous(breaks = seq(0, ceiling(max(n_tf_per_cell$n_tf)), by = 200)) +
        scale_y_continuous(breaks = seq(0, ceiling(max(n_tf_per_cell$mean_umi)), by = 1))
ggMarginal(p, type = "density", col = "brown1", fill = "brown1", alpha = 0.5)

tf_barcodes_f <- tf_barcodes[umi_count > 1]
n_tf_per_cell_f <- tf_barcodes_f[, .(n_tf       = uniqueN(tf_name), 
                                     mean_umi   = as.numeric(mean(umi_count)), 
                                     median_umi = as.numeric(median(umi_count))), by = .(cell_barcode)]

p <- ggplot(n_tf_per_cell_f, aes(n_tf, mean_umi)) +
        geom_point(shape = 19, alpha = 0.4, size = 0.8, col = "royalblue") +
        theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
        theme(axis.title = element_text(size = 20, face = "bold", family = "Arial")) +
        theme(plot.title = element_text(size = 20, face = "bold.italic", family = "Arial")) +
        theme(axis.text = element_text(size = 16, face = "bold")) +
        scale_x_continuous(breaks = replace(seq(0, ceiling(max(n_tf_per_cell_f$n_tf)), by = 200), 1, 10)) +
        scale_y_continuous(breaks = seq(0, ceiling(max(n_tf_per_cell_f$mean_umi)), by = 5)) +
        geom_vline(xintercept = 10, color = "tomato", linetype = "dashed", linewidth = 0.8)
ggMarginal(p, type = "density", col = "brown1", fill = "brown1", alpha = 0.5)

# tf_barcodes_f_plot <- tf_barcodes_f[order(cell_barcode, umi_count), 
#                                     .(tf_rank   = seq_len(.N) / .N, 
#                                       umi_count = umi_count), by = cell_barcode]

bin_values <- c(seq(0, 50, 5), seq(60, 100 , 10), seq(200, 500, 100), Inf)
bin_labels <- c("0~5",     "5~10",    "10~15",   "15~20", 
                "20~25",   "25~30",   "30~35",   "35~40",
                "40~45",   "45~50",   "50~60",   "60~70",
                "70~80",   "80~90",   "90~100",  "100~200",
                "200~300", "300~400", "400~500", "500+")
n_tf_per_cell_f[, tf_bin := cut(n_tf, breaks = bin_values, labels = bin_labels, right=FALSE)]
tf_barcodes_f_plot <- merge(tf_barcodes_f, n_tf_per_cell_f[, .(cell_barcode, tf_bin)], by = "cell_barcode")
cell_counts <- unique(tf_barcodes_f_plot[, .(cell_barcode, tf_bin)])
cell_counts <- cell_counts[, .(num_cells = .N), by = tf_bin]
cell_counts[, tf_bin_label := paste0(tf_bin, " (n=", num_cells, ")")]
cell_counts[, tf_bin_label := factor(paste0(tf_bin, " (n=", num_cells, ")"),
                                     levels = paste0(bin_labels, " (n=", num_cells[match(bin_labels, tf_bin)], ")"))]
tf_barcodes_f_plot <- merge(tf_barcodes_f_plot, cell_counts[, .(tf_bin, tf_bin_label)], by = "tf_bin")

top10_tfs <- tf_barcodes_f_plot[order(cell_barcode, -umi_count),
                                .(top10_min = if(.N >= 10) min(umi_count[1:10]) else min(umi_count)), by = .(cell_barcode, tf_bin, tf_bin_label)]
top5pct <- tf_barcodes_f_plot[, .(top5_umi = quantile(umi_count, 0.95)), by = .(cell_barcode, tf_bin_label)]

tf_barcodes_f_plot1 <- tf_barcodes_f_plot[order(cell_barcode, umi_count),
                                          .(tf_rank   = seq_len(.N) / .N, 
                                            umi_count = umi_count), by = .(cell_barcode, tf_bin_label)]
ggplot(tf_barcodes_f_plot1, aes(x = tf_rank, y = umi_count, group = cell_barcode)) +
    geom_line(alpha = 0.4, color = "royalblue", linewidth = 0.3) +
    geom_hline(data = top10_tfs, aes(yintercept = top10_min, group = cell_barcode),
               color = "red", linetype = "dashed", linewidth = 0.1, alpha = 0.4) +
    facet_wrap(~tf_bin_label, scales = "free", ncol = 5) +
    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
    theme(axis.title = element_text(size = 20, face = "bold", family = "Arial")) +
    theme(plot.title = element_text(size = 20, face = "bold.italic", family = "Arial")) +
    theme(axis.text = element_text(size = 16, face = "bold")) +
    labs(x = "TF rank (per cell)", y = "UMI count")

tf_barcodes_f_plot2 <- tf_barcodes_f_plot[order(cell_barcode, umi_count),
                                          .(tf_rank   = seq_len(.N), 
                                            umi_count = umi_count), by = .(cell_barcode, tf_bin)]
ggplot(tf_barcodes_f_plot2, aes(x = tf_rank, y = umi_count, group = cell_barcode)) +
    geom_line(alpha = 0.4, color = "royalblue", linewidth = 0.3) +
    facet_wrap(~tf_bin, scales = "free", ncol = 5) +
    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
    theme(axis.title = element_text(size = 20, face = "bold", family = "Arial")) +
    theme(plot.title = element_text(size = 20, face = "bold.italic", family = "Arial")) +
    theme(axis.text = element_text(size = 16, face = "bold")) +
    labs(x = "TF rank (per cell)", y = "UMI count")