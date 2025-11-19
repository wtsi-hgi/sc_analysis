#!/usr/bin/env Rscript
quiet_library <- function(pkg) { suppressMessages(suppressWarnings(library(pkg, character.only = TRUE))) }
package_1 <- c("tidyverse", "data.table", "future.apply", "ggplot2", "ggrepel", "optparse", "corrplot", "purrr", "GenomicRanges", "EnsDb.Hsapiens.v86")
package_2 <- c("Seurat", "Signac", "biovizBase", "glmGamPoi", "monocle3", "SeuratWrappers", "clusterProfiler", "org.Hs.eg.db")
packages <- c(package_1, package_2)
invisible(lapply(packages, quiet_library))

option_list <- list(make_option(c("-s", "--sample_sheet"), type = "character", help = "path of the sample sheet",    default = NULL),
                    make_option(c("-b", "--tf_barcode"),   type = "character", help = "path of the TF barcode file", default = NULL),
                    make_option(c("-o", "--output_dir"),   type = "character", help = "output directory",            default = getwd()),
                    make_option(c("-p", "--prefix"),       type = "character", help = "output prefix",               default = NULL),
                    make_option(c("-t", "--threads"),      type = "integer",   help = "number of threads",           default = 1))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if(length(commandArgs(trailingOnly = TRUE)) == 0)
{
    print_help(opt_parser)
    quit(status = 1)
}

# -- check options -- #
if(is.null(opt$sample_sheet))         stop("-s, path of the sample sheet is required!", call. = FALSE)
if(is.null(opt$prefix))               stop("-p, output prefix is required!", call. = FALSE)

#-- function --#
create_overlap_barcode_plot <- function(list_barcodes, plot_file) {
    n <- length(list_barcodes)
    overlap_barcodes <- matrix(0, n, n, dimnames = list(names(list_barcodes), names(list_barcodes)))
    diag(overlap_barcodes) <- 0
    invisible(combn(names(list_barcodes), 2, function(p) {
        i <- p[1]; j <- p[2]
        overlap_barcodes[i, j] <<- length(intersect(list_barcodes[[i]], list_barcodes[[j]]))
        overlap_barcodes[j, i] <<- overlap_barcodes[i, j]
    }))

    png(plot_file, width = 1200, height = 1200, units = "px", res = 200)
    corrplot(overlap_barcodes, method = "color", type = "full", tl.col = "black", addCoef.col = "black", number.cex = 0.5, col = col_palette, cl.pos = "n", is.corr = FALSE)
    dev.off()
}

run_sample_qc <- function(h5file, fragment_file, sample_id) {
    message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "Sample QC: ", sample_id, " --> creating the seurat object ...")
    data <- Read10X_h5(h5file)

    obj <- CreateSeuratObject(counts = data$`Gene Expression`, assay = "RNA")
    obj[["ATAC"]] <- CreateChromatinAssay(counts = data$Peaks, fragments = fragment_file, annotation = annotations, sep = c(":", "-"))
    DefaultAssay(obj) <- "RNA"

    rm(data)
    invisible(gc(verbose = FALSE))

    message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "Sample QC: ", sample_id, " --> QC RNA ...")
    obj$qc_rna_status <- "failed"
    obj$qc_atac_status <- "failed"
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    
    nFeature_low  <- 200
    nFeature_high <- quantile(obj$nFeature_RNA, 0.99)

    nCount_low  <- 1000
    nCount_high <- quantile(obj$nCount_RNA, 0.99)

    qc_rna_obj <- subset(obj, subset = nFeature_RNA > nFeature_low & 
                                       nFeature_RNA < nFeature_high &
                                       nCount_RNA > nCount_low &
                                       nCount_RNA < nCount_high &
                                       percent.mt < 10)
    qc_rna_obj$qc_rna_status <- "passed"
    obj$qc_rna_status[colnames(qc_rna_obj)] <- "passed"

    plot_file <- file.path(qc_dir, paste0(sample_id, ".qc_rna_violin.png"))
    p <- VlnPlot(obj, 
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                 group.by = "qc_rna_status", 
                 cols = c("skyblue","yellowgreen"),
                 pt.size = 0, 
                 ncol = 3)
    ggsave(plot_file, p, width = 10, height = 8)

    message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "Sample QC: ", sample_id, " --> QC ATAC ...")
    DefaultAssay(qc_rna_obj) <- "ATAC"
    qc_rna_obj <- NucleosomeSignal(qc_rna_obj)
    qc_rna_obj <- TSSEnrichment(qc_rna_obj)

    nCount_ATAC_low  <- 1000
    nCount_ATAC_high <- quantile(qc_rna_obj$nCount_ATAC, 0.99)

    tss_enrichment_low  <- 2

    qc_atac_obj <- subset(qc_rna_obj, subset = nCount_ATAC > nCount_ATAC_low & 
                                               nCount_ATAC < nCount_ATAC_high & 
                                               TSS.enrichment > tss_enrichment_low)

    qc_atac_obj$qc_atac_status <- "passed"
    qc_rna_obj$qc_atac_status[colnames(qc_atac_obj)] <- "passed"

    plot_file <- file.path(qc_dir, paste0(sample_id, ".qc_atac_violin.png"))
    p <- VlnPlot(qc_rna_obj, 
                 features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
                 group.by = "qc_atac_status", 
                 cols = c("skyblue","yellowgreen"),
                 pt.size = 0, 
                 ncol = 3)
    ggsave(plot_file, p, width = 10, height = 8)
    
    dt_summary <- data.table(sample_id                  = c(sample_id, sample_id, sample_id),
                             qc_type                    = c("raw", "after_qc_rna", "after_qc_atac"),
                             n_reads_rna                = c(sum(obj@meta.data$nCount_RNA), 
                                                            sum(qc_rna_obj@meta.data$nCount_RNA),
                                                            sum(qc_atac_obj@meta.data$nCount_RNA)),
                             n_reads_atac               = c(sum(obj@meta.data$nCount_ATAC), 
                                                            sum(qc_rna_obj@meta.data$nCount_ATAC),
                                                            sum(qc_atac_obj@meta.data$nCount_ATAC)),
                             n_cells                    = c(nrow(obj@meta.data), 
                                                            nrow(qc_rna_obj@meta.data),
                                                            nrow(qc_atac_obj@meta.data)),
                             n_features_per_cell_median = c(median(obj@meta.data$nFeature_RNA), 
                                                            median(qc_rna_obj@meta.data$nFeature_RNA),
                                                            median(qc_atac_obj@meta.data$nFeature_RNA)),
                             n_reads_per_cell_median    = c(median(obj@meta.data$nCount_RNA), 
                                                            median(qc_rna_obj@meta.data$nCount_RNA),
                                                            median(qc_atac_obj@meta.data$nCount_RNA)),
                             pct_mito_per_cell_median   = c(median(obj@meta.data$percent.mt), 
                                                            median(qc_rna_obj@meta.data$percent.mt),
                                                            median(qc_atac_obj@meta.data$percent.mt)),
                             tss_enrichment_median      = c(NA, 
                                                            median(qc_rna_obj@meta.data$TSS.enrichment),
                                                            median(qc_atac_obj@meta.data$TSS.enrichment)),
                             nucleosome_signal_median   = c(NA,
                                                            median(qc_rna_obj@meta.data$nucleosome_signal),
                                                            median(qc_atac_obj@meta.data$nucleosome_signal)))

    rm(obj)
    rm(qc_rna_obj)
    invisible(gc(verbose = FALSE))

    message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "Sample QC: ", sample_id, " --> normailizing counts ...")
    qc_atac_obj <- RunTFIDF(qc_atac_obj)
    qc_atac_obj <- FindTopFeatures(qc_atac_obj, min.cutoff = "q0")

    DefaultAssay(qc_atac_obj) <- "RNA"
    qc_atac_obj <- SCTransform(qc_atac_obj, verbose = FALSE)

    return(list(qc_seurat_obj = qc_atac_obj, dt_summary = dt_summary))
}


# -- inputs -- #
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "Reading input files ...")

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"

samples <- fread(opt$sample_sheet, header = TRUE, sep = "\t")

col_palette <- colorRampPalette(c("white", "brown1"))(100)

# -- outputs -- #
if(!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)
setwd(opt$output_dir)

sample_prefix <- opt$prefix

qc_dir <- "1_qc_results"
if (!dir.exists(qc_dir)) dir.create(qc_dir)

#-- processing --#
list_barcodes <- list()
list_h5_paths <- list()
list_fragments <- list()
for(i in seq_along(samples$sample_id)) {
    list_barcodes[[samples$sample_id[i]]] <- fread(paste0(samples$cellranger_dir[i], "/per_barcode_metrics.csv"), header = TRUE, sep = ",", select = c("barcode", "is_cell"))
    list_barcodes[[samples$sample_id[i]]] <- list_barcodes[[samples$sample_id[i]]][is_cell == 1]$barcode

    list_h5_paths[[samples$sample_id[i]]] <- paste0(samples$cellranger_dir[i], "/filtered_feature_bc_matrix.h5")
    list_fragments[[samples$sample_id[i]]] <- paste0(samples$cellranger_dir[i], "/atac_fragments.tsv.gz")
}

# 1. barcode checking
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "Checking barcodes ...")
plot_file <- file.path(qc_dir, paste0(sample_prefix, ".overlapped_barcodes.png"))
create_overlap_barcode_plot(list_barcodes, plot_file)

rm(list_barcodes)
invisible(gc(verbose = FALSE))

# 2. QC for each sample
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "Sample QC ...")

# parallel processing needs lots of MEM, like ~600GB here, not effeicient for the job, so using loop here 
# plan(multicore, workers = opt$threads) 
# future_res <- future_mapply(run_sample_qc, h5file = list_h5_paths, fragment_file = list_fragments, sample_id = names(list_h5_paths), SIMPLIFY = FALSE)
# qc_seurat_objs <- map(future_res, "qc_seurat_obj")
# qc_summary <- rbindlist(map(future_res, "dt_summary"))

qc_seurat_objs <- list()
qc_summary <- data.table()
for(i in seq_along(samples$sample_id)) {
    idx <- samples$sample_id[i]
    qc_results <- run_sample_qc(list_h5_paths[[idx]], list_fragments[[idx]], idx)
    qc_seurat_objs[[idx]] <- qc_results$qc_seurat_obj
    qc_summary <- rbindlist(list(qc_summary, qc_results$dt_summary), use.names = TRUE, fill = TRUE)
}

fwrite(qc_summary, file = file.path(qc_dir, paste0(sample_prefix, ".qc_summary.tsv")), sep = "\t")

list_qc_barcodes <- lapply(qc_seurat_objs, colnames)
plot_file <- file.path(qc_dir, paste0(sample_prefix, ".qc_overlapped_barcodes.png"))
create_overlap_barcode_plot(list_qc_barcodes, plot_file)

dt_cell_barcodes <- as.data.table(stack(list_qc_barcodes))
dt_cell_barcodes <- dt_cell_barcodes[, .(count = .N, samples = paste(ind, collapse = ",")), by = values]
dt_cell_barcodes_overlap <- dt_cell_barcodes[count > 1]
setorder(dt_cell_barcodes_overlap, -count, values)
fwrite(dt_cell_barcodes_overlap, file = file.path(qc_dir, paste0(sample_prefix, ".qc_overlapped_barcodes.tsv")), sep = "\t")

qc_seurat_objs <- lapply(qc_seurat_objs, function(obj) {
    keep <- !colnames(obj) %in% dt_cell_barcodes_overlap$values
    obj[, keep]
})
saveRDS(qc_seurat_objs, file = file.path(qc_dir, paste0(sample_prefix, ".qc_seurat_objs.unique_cells.rds")))

rm(list_qc_barcodes)
invisible(gc(verbose = FALSE))

# 3. integration
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "Removing batch effects ...")
features <- SelectIntegrationFeatures(object.list = qc_seurat_objs, nfeatures = 3000)
qc_seurat_objs <- PrepSCTIntegration(object.list = qc_seurat_objs, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = qc_seurat_objs, normalization.method = "SCT", anchor.features = features)
integrated_obj <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

rm(qc_seurat_objs)
invisible(gc(verbose = FALSE))

old_names <- colnames(integrated_obj)
new_names <- sub("-1_.*$", "", old_names)
dup_cells <- duplicated(new_names)
integrated_obj <- integrated_obj[, !dup_cells]
colnames(integrated_obj) <- new_names[!dup_cells]
cells_in_obj <- colnames(integrated_obj)

tf_barcodes_umi_top10 <- fread(opt$tf_barcode, sep = "\t", header = T)
tf_barcodes_umi_top10 <- tf_barcodes_umi_top10[cell_barcode %in% cells_in_obj]
missing_cells <- setdiff(cells_in_obj, tf_barcodes_umi_top10$cell_barcode)

integrated_obj <- subset(integrated_obj, cells = setdiff(colnames(integrated_obj), missing_cells))

tf_mat <- dcast(tf_barcodes_umi_top10, tf_name ~ cell_barcode, value.var = "umi_count", fill = 0)
tf_mat <- tf_mat[!is.na(tf_name)]
tf_mat <- as.data.frame(tf_mat)
rownames(tf_mat) <- tf_mat$tf_name
tf_mat$tf_name <- NULL
tf_mat <- as.matrix(tf_mat[, colnames(integrated_obj)])

tf_assay <- CreateAssayObject(data = tf_mat)
integrated_obj[["TF"]] <- tf_assay
integrated_obj[["TF"]]@counts <- tf_mat
DefaultAssay(integrated_obj) <- "TF"
integrated_obj <- SCTransform(integrated_obj, assay = "TF", new.assay.name = "SCT_TF", verbose = FALSE)

DefaultAssay(integrated_obj) <- "integrated"
integrated_obj <- FindVariableFeatures(integrated_obj, assay = "integrated", selection.method = "vst", nfeatures = 3000)
integrated_obj <- ScaleData(integrated_obj, assay = "integrated", verbose = FALSE)
integrated_obj <- RunPCA(integrated_obj, assay = "integrated", npcs = 30, verbose = FALSE)
integrated_obj <- RunUMAP(integrated_obj, assay = "integrated", dims = 1:30)
integrated_obj <- FindNeighbors(integrated_obj, assay = "integrated", dims = 1:30)
integrated_obj <- FindClusters(integrated_obj, algorithm = 1, assay = "integrated", resolution = 0.5)

saveRDS(integrated_obj, file = file.path(qc_dir, paste0(sample_prefix, ".integrated_obj.rds")))
DimPlot(integrated_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

FeaturePlot(integrated_obj, features = "PC_1")

# 4. TF analysis
clusters_ident <- levels(integrated_obj@meta.data$seurat_clusters)
tf_expressions <- list()
for (cl in clusters_ident) {
    cells_in_cluster <- WhichCells(integrated_obj, idents = cl)
    tf_data <- GetAssayData(integrated_obj, assay = "SCT_TF", slot = "data")[, cells_in_cluster]
    
    pct_expressed <- rowSums(tf_data > 0) / length(cells_in_cluster)
    expressed_tfs <- names(pct_expressed[pct_expressed > 0.01])
  
    mean_vals <- sapply(expressed_tfs, function(tf) { 
        vals <- tf_data[tf, ]
        vals <- vals[vals > 0]
        if (length(vals) == 0) NA else mean(vals)
    })
  
    dt <- data.table(tf         = expressed_tfs,
                     mean_val   = mean_vals[expressed_tfs],
                     exp_ratio  = pct_expressed[expressed_tfs],
                     cluster    = cl,
                     n_cells    = length(cells_in_cluster))
    setorder(dt, -exp_ratio)
    tf_expressions[[cl]] <- dt
}

for (i in 1:length(clusters_ident)) {
    dt <- tf_expressions[[i]]
    top_tfs <- dt[order(-exp_ratio)][1:min(10, .N)]
    
    plot_file <- paste0("morf10_tf.integration_cl", unique(dt$cluster), ".tf_profile.png")

    p <- ggplot(dt, aes(x = mean_val, y = exp_ratio)) +
                geom_point(shape = 16, color = "lightgrey") +
                geom_point(data = top_tfs, color = "red", size = 3) +
                geom_segment(data = top_tfs, aes(x = mean_val, y = exp_ratio, xend = mean_val, yend = exp_ratio), alpha = 0) +
                geom_text_repel(data = top_tfs, aes(label = tf), size = 3.5, point.padding = 0.2, box.padding = 0.5, 
                                segment.size = 0.4, segment.color = "black", min.segment.length = 0,
                                arrow = arrow(length = unit(0.2, "cm"), type = "open", ends = "last")) +
                theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                theme(axis.title = element_text(size = 10, face = "bold", family = "Arial")) +
                theme(plot.title = element_text(size = 10, face = "bold.italic", family = "Arial")) +
                theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                labs(x = "Mean Normalised UMI", y = "Ratio of Expressed Cells", 
                    title = paste0("Cluster ", unique(dt$cluster), " (n=", unique(dt$n_cells), ")"))

    png(plot_file, width = 1000, height = 1000, units = "px", res = 200)
    print(p)
    dev.off()
}

# 5. marker genes and TFs for each cluster
gene_marks_clusters <- lapply(clusters_ident, function(cl) {
    FindMarkers(integrated_obj,
    ident.1 = cl,
    assay = "SCT",
    only.pos = TRUE,
    logfc.threshold = 0.25,
    min.pct = 0.1)
})
names(gene_marks_clusters) <- paste0("cluster_", clusters_ident)
sapply(gene_marks_clusters, nrow)
gene_marks_clusters <- lapply(gene_marks_clusters, function(df) {if (!is.null(df)) as.data.table(df, keep.rownames = "gene") else NULL})

gene_marks_clusters_sig <- lapply(gene_marks_clusters, function(x) {if (!is.null(x)) x[p_val <= 0.05] else NULL})
sapply(gene_marks_clusters_sig, nrow)

dt_gene_marks_clusters_sig <- do.call(rbind, lapply(names(gene_marks_clusters_sig), function(cl) {
    dt <- gene_marks_clusters_sig[[cl]]
    dt$cluster <- cl
    return(dt)
}))

fwrite(dt_gene_marks_clusters_sig, "morf10_tf.integration.cluster_sig_genes.tsv", quote = F, sep = "\t")

tf_marks_clusters <- lapply(clusters_ident, function(cl) {
  FindMarkers(integrated_obj,
              ident.1 = cl,
              assay = "SCT_TF",
              only.pos = TRUE,
              logfc.threshold = 0.25,
              min.pct = 0.01)
})
names(tf_marks_clusters) <- paste0("cluster_", clusters_ident)
sapply(tf_marks_clusters, nrow)
tf_marks_clusters <- lapply(tf_marks_clusters, function(df) {if (!is.null(df)) as.data.table(df, keep.rownames = "tf") else NULL})

tf_marks_clusters_sig <- lapply(tf_marks_clusters, function(x) {if (!is.null(x)) x[p_val <= 0.05] else NULL})
sapply(tf_marks_clusters_sig, nrow)

dt_tf_marks_clusters_sig <- do.call(rbind, lapply(names(tf_marks_clusters_sig), function(cl) {
    dt <- tf_marks_clusters_sig[[cl]]
    dt$cluster <- cl
    return(dt)
}))

fwrite(dt_tf_marks_clusters_sig, "morf10_tf.integration.cluster_sig_tfs.tsv", quote = F, sep = "\t")

for (i in 1:length(clusters_ident)) {
    dt <- tf_expressions[[i]]
    top_tfs <- dt[tf %in% tf_marks_clusters_sig[[i]]$tf]
    
    plot_file <- paste0("morf10_tf.integration_cl", unique(dt$cluster), ".tf_profile_sig.png")

    p <- ggplot(dt, aes(x = mean_val, y = exp_ratio)) +
                geom_point(shape = 16, color = "lightgrey") +
                geom_point(data = top_tfs, color = "red", size = 3) +
                geom_segment(data = top_tfs, aes(x = mean_val, y = exp_ratio, xend = mean_val, yend = exp_ratio), alpha = 0) +
                geom_text_repel(data = top_tfs, aes(label = tf), size = 3.5, point.padding = 0.2, box.padding = 0.5, 
                                segment.size = 0.4, segment.color = "black", min.segment.length = 0,
                                arrow = arrow(length = unit(0.2, "cm"), type = "open", ends = "last")) +
                theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                theme(axis.title = element_text(size = 10, face = "bold", family = "Arial")) +
                theme(plot.title = element_text(size = 10, face = "bold.italic", family = "Arial")) +
                theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                labs(x = "Mean Normalised UMI", y = "Ratio of Expressed Cells", 
                    title = paste0("Cluster ", unique(dt$cluster), " (n=", unique(dt$n_cells), ")"))

    png(plot_file, width = 1000, height = 1000, units = "px", res = 200)
    print(p)
    dev.off()
}









gene_marks_cluster17 <- FindMarkers(integrated_obj,
                                    ident.1 = 17,
                                    assay = "SCT",       
                                    only.pos = TRUE,        
                                    logfc.threshold = 0.25,
                                    min.pct = 0.10)
tf_marks_cluster17 <- FindMarkers(integrated_obj,
                                  ident.1 = 17,
                                  assay = "SCT_TF",       
                                  only.pos = TRUE,        
                                  logfc.threshold = 0.25,
                                  min.pct = 0.01)

gene_marks_cluster17_sig <- gene_marks_cluster17 %>%
                            rownames_to_column(var = "gene") %>%
                            dplyr::filter(p_val_adj < 0.05 & avg_log2FC > 0.25)

genes_entrez <- bitr(gene_marks_cluster17_sig$gene, fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)
ego <- enrichGO(gene          = genes_entrez$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff  = 0.05,
                readable      = TRUE)

barplot(ego, showCategory = 20)
dotplot(ego, showCategory = 20)
emapplot(ego)

go_results <- as.data.frame(ego)









DefaultAssay(integrated_obj) <- "SCT"
cds <- as.cell_data_set(integrated_obj)
reducedDim(cds, "UMAP") <- integrated_obj@reductions$umap@cell.embeddings
cds@clusters$UMAP$clusters <- as.factor(integrated_obj$seurat_clusters)
names(cds@clusters$UMAP$clusters) <- colnames(cds)
cds@clusters$UMAP$partitions <- rep(1, ncol(cds))
names(cds@clusters$UMAP$partitions) <- colnames(cds)
cds <- learn_graph(cds)

root_cells <- colnames(cds)[cds@clusters$UMAP$clusters == "17"]
cds <- order_cells(cds, root_cells = root_cells)
plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = F)

save.image("morf10_tf.integration.RData")

# integrated_obj <- AddMetaData(integrated_obj, metadata = pseudotime(cds)[colnames(integrated_obj)], col.name = "pseudotime")

DefaultAssay(integrated_obj) <- "integrated"
hvg_genes <- VariableFeatures(integrated_obj)

cds <- detect_genes(cds)
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 100))

genes_to_test <- intersect(expressed_genes, hvg_genes)
cds_sub <- cds[genes_to_test, ]

cds_fit <- fit_models(cds_sub, model_formula_str = "~pseudotime", cores = 1, clean_model = TRUE)
fit_coefs <- coefficient_table(cds_fit)
pseudo_res <- subset(fit_coefs, term == "pseudotime")
pseudo_res$direction <- ifelse(pseudo_res$estimate > 0, "up", "down")
pseudo_res$score <- pseudo_res$estimate * -log10(pseudo_res$p_value)
top_up <- head(pseudo_res[pseudo_res$direction == "up", ][order(-pseudo_res$score), ], 1000)
top_down <- head(pseudo_res[pseudo_res$direction == "down", ][order(-pseudo_res$score), ], 1000)





