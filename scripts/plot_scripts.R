library(data.table)
library(ggplot2)
library(ggExtra)
library(reshape2)
library(scales)

input_dir <- "/lustre/scratch127/gengen/projects_v2/dual_tf/MORF_10/Multiome/analysis_by_hgi/1_qc_results"

qc_res <- fread(file.path(input_dir, "morf10.qc_summary.tsv"))
qc_melted <- as.data.table(melt(qc_res,
                  id.vars = c("sample_id", "qc_type"),
                  measure.vars = c("n_reads_rna", "n_reads_atac", "n_cells", "n_features_per_cell_median"),
                  variable.name = "metric",
                  value.name = "value"))
qc_melted[, qc_type := factor(qc_type, levels = c("raw", "after_qc_rna", "after_qc_atac"))]

ggplot(qc_melted, aes(x = sample_id, y = value, fill = qc_type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +
    facet_wrap(~metric, scales = "free_y") +
    scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
    labs(x = "Sample", y = "Value", fill = "QC Type") +
    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
    theme(strip.text = element_text(size = 16, face = "bold")) +
    theme(axis.title = element_text(size = 20, face = "bold", family = "Arial")) +
    theme(plot.title = element_text(size = 20, face = "bold.italic", family = "Arial")) +
    theme(axis.text = element_text(size = 16, face = "bold")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

input_dir <- "/lustre/scratch127/gengen/projects_v2/dual_tf/MORF_10/Multiome/analysis_by_hgi/2_tf_barcodes"

# 1. filter tf barcode with raw read count < 2
# 2. umi should be assigned to only one tf, so keep the tf with highest umis
# 3. filter tf barcode with umi count < 2
tf_barcodes <- fread(file.path(input_dir, "morf10_tf.barcode_out.tsv.gz"))
tf_count_a1 <- tf_barcodes[count > 1]
tf_count_a1 <- tf_count_a1[order(-count), .SD[1], by = .(cell_barcode, read_umi)]
tf_count_a1_umi <- tf_count_a1[, .(umi_count = uniqueN(read_umi)), by = .(cell_barcode, tf_name)]
tf_count_a1_umi_a1 <- tf_count_a1_umi[umi_count > 1]

fwrite(tf_count_a1_umi_a1, "morf10_tf.umi_filtered.tsv", row.names = F, sep = "\t")

tf_count_a1_umi_n_tf <- tf_count_a1_umi[, .(n_tf       = uniqueN(tf_name), 
                                            mean_umi   = as.numeric(mean(umi_count)), 
                                            median_umi = as.numeric(median(umi_count))), by = .(cell_barcode)]
p <- ggplot(tf_count_a1_umi_n_tf, aes(n_tf, mean_umi)) +
        geom_point(shape = 19, alpha = 0.4, size = 0.8, col = "royalblue") +
        theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
        theme(axis.title = element_text(size = 20, face = "bold", family = "Arial")) +
        theme(plot.title = element_text(size = 20, face = "bold.italic", family = "Arial")) +
        theme(axis.text = element_text(size = 16, face = "bold")) +
        scale_x_continuous(breaks = replace(seq(0, ceiling(max(tf_count_a1_umi_n_tf$n_tf)), by = 30), 1, 10)) +
        scale_y_continuous(breaks = seq(0, ceiling(max(tf_count_a1_umi_n_tf$mean_umi)), by = 5)) +
        geom_vline(xintercept = 10, color = "tomato", linetype = "dashed", linewidth = 0.8)
ggMarginal(p, type = "density", col = "brown1", fill = "brown1", alpha = 0.5)

tf_count_a1_umi_a1_n_tf <- tf_count_a1_umi_a1[, .(n_tf       = uniqueN(tf_name), 
                                                  mean_umi   = as.numeric(mean(umi_count)), 
                                                  median_umi = as.numeric(median(umi_count))), by = .(cell_barcode)]
p <- ggplot(tf_count_a1_umi_a1_n_tf, aes(n_tf, mean_umi)) +
        geom_point(shape = 19, alpha = 0.4, size = 0.8, col = "royalblue") +
        theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
        theme(axis.title = element_text(size = 20, face = "bold", family = "Arial")) +
        theme(plot.title = element_text(size = 20, face = "bold.italic", family = "Arial")) +
        theme(axis.text = element_text(size = 16, face = "bold")) +
        scale_x_continuous(breaks = replace(seq(0, ceiling(max(tf_count_a1_umi_a1_n_tf$n_tf)), by = 30), 1, 10)) +
        scale_y_continuous(breaks = seq(0, ceiling(max(tf_count_a1_umi_a1_n_tf$mean_umi)), by = 10)) +
        geom_vline(xintercept = 10, color = "tomato", linetype = "dashed", linewidth = 0.8)
ggMarginal(p, type = "density", col = "brown1", fill = "brown1", alpha = 0.5)

bin_values <- c(seq(0, 50, 5), seq(60, 100 , 10), Inf)
bin_labels <- c("0~5",     "5~10",    "10~15",   "15~20", 
                "20~25",   "25~30",   "30~35",   "35~40",
                "40~45",   "45~50",   "50~60",   "60~70",
                "70~80",   "80~90",   "90~100",  "100+")
tf_count_a1_umi_n_tf[, tf_bin := cut(n_tf, breaks = bin_values, labels = bin_labels, right=FALSE)]
tf_count_a1_umi_plot <- merge(tf_count_a1_umi, tf_count_a1_umi_n_tf[, .(cell_barcode, tf_bin)], by = "cell_barcode")
cell_counts <- unique(tf_count_a1_umi_plot[, .(cell_barcode, tf_bin)])
cell_counts <- cell_counts[, .(num_cells = .N), by = tf_bin]
cell_counts[, tf_bin_label := paste0(tf_bin, " (n=", num_cells, ")")]
cell_counts[, tf_bin_label := factor(paste0(tf_bin, " (n=", num_cells, ")"),
                                     levels = paste0(bin_labels, " (n=", num_cells[match(bin_labels, tf_bin)], ")"))]
tf_count_a1_umi_plot <- merge(tf_count_a1_umi_plot, cell_counts[, .(tf_bin, tf_bin_label)], by = "tf_bin")

top10_tfs <- tf_count_a1_umi_plot[order(cell_barcode, -umi_count),
                                  .(top10_min = if(.N >= 10) min(umi_count[1:10]) else min(umi_count)), by = .(cell_barcode, tf_bin, tf_bin_label)]
#top5pct <- tf_barcodes_f_plot[, .(top5_umi = quantile(umi_count, 0.95)), by = .(cell_barcode, tf_bin_label)]

tf_count_a1_umi_plot1 <- tf_count_a1_umi_plot[order(cell_barcode, umi_count),
                                              .(tf_rank   = seq_len(.N) / .N, 
                                                umi_count = umi_count), by = .(cell_barcode, tf_bin_label)]
ggplot(tf_count_a1_umi_plot1, aes(x = tf_rank, y = umi_count, group = cell_barcode)) +
    geom_line(alpha = 0.4, color = "royalblue", linewidth = 0.3) +
    geom_hline(data = top10_tfs, aes(yintercept = top10_min, group = cell_barcode),
               color = "red", linetype = "dashed", linewidth = 0.1, alpha = 0.4) +
    facet_wrap(~tf_bin_label, scales = "free", ncol = 4) +
    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
    theme(axis.title = element_text(size = 10, face = "bold", family = "Arial")) +
    theme(plot.title = element_text(size = 10, face = "bold.italic", family = "Arial")) +
    theme(axis.text = element_text(size = 8, face = "bold")) +
    labs(x = "TF rank (per cell)", y = "UMI count")

tf_count_a1_umi_plot2 <- tf_count_a1_umi_plot[order(cell_barcode, umi_count),
                                              .(tf_rank   = seq_len(.N), 
                                                umi_count = umi_count), by = .(cell_barcode, tf_bin_label)]
ggplot(tf_count_a1_umi_plot2, aes(x = tf_rank, y = umi_count, group = cell_barcode)) +
    geom_line(alpha = 0.4, color = "royalblue", linewidth = 0.3) +
    facet_wrap(~tf_bin_label, scales = "free", ncol = 4) +
    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
    theme(axis.title = element_text(size = 10, face = "bold", family = "Arial")) +
    theme(plot.title = element_text(size = 10, face = "bold.italic", family = "Arial")) +
    theme(axis.text = element_text(size = 8, face = "bold")) +
    labs(x = "TF rank (per cell)", y = "UMI count")

tf_count_a1_umi_a1_n_tf[, tf_bin := cut(n_tf, breaks = bin_values, labels = bin_labels, right=FALSE)]
tf_count_a1_umi_a1_plot <- merge(tf_count_a1_umi_a1, tf_count_a1_umi_a1_n_tf[, .(cell_barcode, tf_bin)], by = "cell_barcode")
cell_counts <- unique(tf_count_a1_umi_a1_plot[, .(cell_barcode, tf_bin)])
cell_counts <- cell_counts[, .(num_cells = .N), by = tf_bin]
cell_counts[, tf_bin_label := paste0(tf_bin, " (n=", num_cells, ")")]
cell_counts[, tf_bin_label := factor(paste0(tf_bin, " (n=", num_cells, ")"),
                                     levels = paste0(bin_labels, " (n=", num_cells[match(bin_labels, tf_bin)], ")"))]
tf_count_a1_umi_a1_plot <- merge(tf_count_a1_umi_a1_plot, cell_counts[, .(tf_bin, tf_bin_label)], by = "tf_bin")

top10_tfs <- tf_count_a1_umi_a1_plot[order(cell_barcode, -umi_count),
                                    .(top10_min = if(.N >= 10) min(umi_count[1:10]) else min(umi_count)), by = .(cell_barcode, tf_bin, tf_bin_label)]

tf_count_a1_umi_a1_plot1 <- tf_count_a1_umi_a1_plot[order(cell_barcode, umi_count),
                                                    .(tf_rank   = seq_len(.N) / .N, 
                                                      umi_count = umi_count), by = .(cell_barcode, tf_bin_label)]
ggplot(tf_count_a1_umi_a1_plot1, aes(x = tf_rank, y = umi_count, group = cell_barcode)) +
    geom_line(alpha = 0.4, color = "royalblue", linewidth = 0.3) +
    geom_hline(data = top10_tfs, aes(yintercept = top10_min, group = cell_barcode),
               color = "red", linetype = "dashed", linewidth = 0.1, alpha = 0.4) +
    facet_wrap(~tf_bin_label, scales = "free", ncol = 4) +
    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
    theme(axis.title = element_text(size = 10, face = "bold", family = "Arial")) +
    theme(plot.title = element_text(size = 10, face = "bold.italic", family = "Arial")) +
    theme(axis.text = element_text(size = 8, face = "bold")) +
    labs(x = "TF rank (per cell)", y = "UMI count")

tf_count_a1_umi_a1_plot2 <- tf_count_a1_umi_a1_plot[order(cell_barcode, umi_count),
                                                    .(tf_rank   = seq_len(.N), 
                                                      umi_count = umi_count), by = .(cell_barcode, tf_bin_label)]
ggplot(tf_count_a1_umi_a1_plot2, aes(x = tf_rank, y = umi_count, group = cell_barcode)) +
    geom_line(alpha = 0.4, color = "royalblue", linewidth = 0.3) +
    facet_wrap(~tf_bin_label, scales = "free", ncol = 4) +
    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
    theme(axis.title = element_text(size = 10, face = "bold", family = "Arial")) +
    theme(plot.title = element_text(size = 10, face = "bold.italic", family = "Arial")) +
    theme(axis.text = element_text(size = 8, face = "bold")) +
    labs(x = "TF rank (per cell)", y = "UMI count")