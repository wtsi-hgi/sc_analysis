#!/usr/bin/env Rscript
quiet_library <- function(pkg) { suppressMessages(suppressWarnings(library(pkg, character.only = TRUE))) }
package_1 <- c("optparse", "tidyverse", "data.table", "ggplot2", "ggrepel")
package_2 <- c("Seurat", "Signac", "SeuratWrappers", "clusterProfiler", "GenomicRanges", "EnsDb.Hsapiens.v86")
packages <- c(package_1, package_2)
invisible(lapply(packages, quiet_library))

option_list <- list(make_option(c("-s", "--sample_id"),  type = "character", help = "sample ID",               default = NULL),
                    make_option(c("-g", "--gex_file"),   type = "character", help = "gene expression h5 file", default = NULL),
                    make_option(c("-a", "--atac_file"),  type = "character", help = "ATAC tsv file",           default = NULL),
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
if(is.null(opt$sample_id)) stop("-s, sample ID is required!", call. = FALSE)
if(is.null(opt$gex_file))  stop("-g, gene expression h5 file is required!", call. = FALSE)
if(is.null(opt$atac_file)) stop("-a, ATAC tsv file is required!", call. = FALSE)

# -- inputs -- #
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"

# -- outputs -- #
if(!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)
setwd(opt$output_dir)

sample_prefix <- ifelse(is.null(opt$prefix), opt$sample_id, opt$prefix)

#-- processing --#
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "Sample QC: ", opt$sample_id, " --> creating the seurat object ...")

data_gex <- Read10X_h5(opt$gex_file)
obj <- CreateSeuratObject(counts = data_gex$`Gene Expression`, assay = "RNA")
obj[["ATAC"]] <- CreateChromatinAssay(counts = data_atac$Peaks, fragments = opt$atac_file, annotation = annotations, sep = c(":", "-"))

rm(data_gex)
invisible(gc(verbose = FALSE))

message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "Sample QC: ", opt$sample_id, " --> QC RNA ...")
DefaultAssay(obj) <- "RNA"
obj$qc_rna_status <- "failed"

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

plot_file <- paste0(sample_prefix, ".qc_rna_violin.png")
p <- VlnPlot(obj, 
             features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
             group.by = "qc_rna_status", 
             cols = c("skyblue","yellowgreen"),
             pt.size = 0, 
             ncol = 3)
ggsave(plot_file, p, width = 10, height = 8)

message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "Sample QC: ", opt$sample_id, " --> QC ATAC ...")
DefaultAssay(qc_rna_obj) <- "ATAC"
obj$qc_atac_status <- "failed"

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

plot_file <- paste0(sample_prefix, ".qc_atac_violin.png")
p <- VlnPlot(qc_rna_obj, 
             features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
             group.by = "qc_atac_status", 
             cols = c("skyblue","yellowgreen"),
             pt.size = 0, 
             ncol = 3)
ggsave(plot_file, p, width = 10, height = 8)

dt_summary <- data.table(sample_id                  = c(opt$sample_id, opt$sample_id, opt$sample_id),
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

message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "Sample QC: ", opt$sample_id, " --> normailizing counts ...")
qc_atac_obj <- RunTFIDF(qc_atac_obj)
qc_atac_obj <- FindTopFeatures(qc_atac_obj, min.cutoff = "q0")

DefaultAssay(qc_atac_obj) <- "RNA"
qc_atac_obj <- SCTransform(qc_atac_obj, verbose = FALSE)

message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "Sample QC: ", opt$sample_id, " --> creating output files ...")
saveRDS(qc_atac_obj, file = paste0(sample_prefix ".qc_obj.rds"))
fwrite(dt_summary, file = paste0(sample_prefix, ".qc_summary.tsv"), sep = "\t")
fwrite(colnames(qc_atac_obj), file = paste0(sample_prefix, ".qc_cell_barcodes.tsv"), sep = "\t", col.names = FALSE)
