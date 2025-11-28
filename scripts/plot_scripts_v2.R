library(data.table)
library(ggplot2)
library(ggExtra)
library(reshape2)
library(scales)

input_dir <- "/lustre/scratch127/gengen/projects_v2/dual_tf/MORF_10/Multiome/analysis_by_hgi/1_qc_scrna"

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

input_dir <- "/lustre/scratch127/gengen/projects_v2/dual_tf/MORF_10/Multiome/analysis_by_hgi/2_qc_sctf"

tf_counts <- fread(file.path(input_dir, "50608.barcode_out.tsv.gz"))
tf_counts <- tf_counts[, .(count = sum(count)), by = .(cell_barcode, tf_name)]
tf_counts <- tf_counts[count > 1]
tf_counts[, cluster := {
    if (.N > 10) {
        ck <- Ckmeans.1d.dp(log2(count), k = 2, y = 1)$cluster
        ck
    } else {
        rep(2L, .N)
    }
}, by = cell_barcode]

tf_counts_f <- tf_counts[cluster == 2]
tf_counts_f[, cluster := NULL]

fwrite(tf_counts_f, "50608.tf_per_cell_kmeans.tsv", row.names = F, sep = "\t")

tf_counts_f_n_tf <- tf_counts_f[, .(n_tf       = uniqueN(tf_name), 
                                    mean_count   = as.numeric(mean(count)), 
                                    median_count = as.numeric(median(count))), by = .(cell_barcode)]
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

bin_values <- c(seq(0, 50, 5), seq(60, 100 , 10), Inf)
bin_labels <- c("0~5",     "5~10",    "10~15",   "15~20", 
                "20~25",   "25~30",   "30~35",   "35~40",
                "40~45",   "45~50",   "50~60",   "60~70",
                "70~80",   "80~90",   "90~100",  "100+")
tf_counts_f_n_tf[, tf_bin := cut(n_tf, breaks = bin_values, labels = bin_labels, right=FALSE)]
tf_counts_f_plot <- merge(tf_counts_f, tf_counts_f_n_tf[, .(cell_barcode, tf_bin)], by = "cell_barcode")
cell_counts <- unique(tf_counts_f_plot[, .(cell_barcode, tf_bin)])
cell_counts <- cell_counts[, .(num_cells = .N), by = tf_bin]
cell_counts[, tf_bin_label := paste0(tf_bin, " (n=", num_cells, ")")]
cell_counts[, tf_bin_label := factor(paste0(tf_bin, " (n=", num_cells, ")"),
                                     levels = paste0(bin_labels, " (n=", num_cells[match(bin_labels, tf_bin)], ")"))]
tf_counts_f_plot <- merge(tf_counts_f_plot, cell_counts[, .(tf_bin, tf_bin_label)], by = "tf_bin")

top10_tfs <- tf_counts_f_plot[order(cell_barcode, -count),
                                  .(top10_min = if(.N >= 10) min(count[1:10]) else min(count)), by = .(cell_barcode, tf_bin, tf_bin_label)]

tf_counts_f_plot1 <- tf_counts_f_plot[order(cell_barcode, count),
                                      .(tf_rank   = seq_len(.N) / .N, 
                                        count = count), by = .(cell_barcode, tf_bin_label)]
ggplot(tf_counts_f_plot1, aes(x = tf_rank, y = count, group = cell_barcode)) +
    geom_line(alpha = 0.4, color = "royalblue", linewidth = 0.3) +
    geom_hline(data = top10_tfs, aes(yintercept = top10_min, group = cell_barcode),
               color = "red", linetype = "dashed", linewidth = 0.1, alpha = 0.4) +
    facet_wrap(~tf_bin_label, scales = "free", ncol = 4) +
    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
    theme(axis.title = element_text(size = 10, face = "bold", family = "Arial")) +
    theme(plot.title = element_text(size = 10, face = "bold.italic", family = "Arial")) +
    theme(axis.text = element_text(size = 8, face = "bold")) +
    labs(x = "TF rank (per cell)", y = "Count")

tf_counts_f_plot2 <- tf_counts_f_plot[order(cell_barcode, count),
                                      .(tf_rank   = seq_len(.N), 
                                        count = count), by = .(cell_barcode, tf_bin_label)]
ggplot(tf_counts_f_plot2, aes(x = tf_rank, y = count, group = cell_barcode)) +
    geom_line(alpha = 0.4, color = "royalblue", linewidth = 0.3) +
    facet_wrap(~tf_bin_label, scales = "free", ncol = 4) +
    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
    theme(axis.title = element_text(size = 10, face = "bold", family = "Arial")) +
    theme(plot.title = element_text(size = 10, face = "bold.italic", family = "Arial")) +
    theme(axis.text = element_text(size = 8, face = "bold")) +
    labs(x = "TF rank (per cell)", y = "Count")

tf_counts_f_boxplot <- melt(tf_counts_f_n_tf, 
                            id.vars = c("cell_barcode", "tf_bin"), 
                            measure.vars = c("mean_count", "median_count"),
                            variable.name = "count_type", value.name = "count_value")
ggplot(tf_counts_f_boxplot, aes(x = tf_bin, y = count_value, fill = count_type)) +
    geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.8), outlier.colour = "darkgrey", linewidth = 0.1) +
    scale_y_log10() +
    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
    theme(axis.title = element_text(size = 10, face = "bold", family = "Arial")) +
    theme(plot.title = element_text(size = 10, face = "bold.italic", family = "Arial")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "TF Bin", y = "Count", title = "Count by TF Bin")

