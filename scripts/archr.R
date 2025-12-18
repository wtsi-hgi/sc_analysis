library(ArchR)
addArchRGenome("hg38")
addArchRThreads(16)
addArchRLocking(locking = TRUE)
set.seed(1)

setwd("/lustre/scratch127/gengen/projects_v2/dual_tf/MORF_10/Multiome/analysis_by_hgi/4_archr")

library(DropletUtils)
sample_id <- "morf10"
i <- 0
for (layer in Layers(integrated_obj[["RNA"]])) {
    message("Exporting ", layer)

    old_names <- colnames(integrated_obj)
    new_names <- paste0(new_names, "-1")
    colnames(integrated_obj) <- new_names

    mat <- GetAssayData(integrated_obj, assay = "RNA", layer = layer)
    mat@x <- as.integer(round(mat@x))

    data    <- mat@x
    indices <- mat@i
    indptr  <- mat@p
    shape   <- as.integer(dim(mat))
    barcodes <- colnames(mat)
    features <- rownames(mat) 

    i <- i + 1
    h5file <- paste0(sample_id, "_b", i, "_qced.h5")

    h5createFile(h5file)
    h5createGroup(h5file, "matrix")
    h5createGroup(h5file, "matrix/features")

    h5write(data,    h5file, "matrix/data")
    h5write(indices, h5file, "matrix/indices")
    h5write(indptr,  h5file, "matrix/indptr")
    h5write(shape,   h5file, "matrix/shape")

    h5write(as.character(colnames(mat)), h5file, "matrix/barcodes")
    h5write(as.character(features), h5file, "matrix/features/id")
    h5write(as.character(features), h5file, "matrix/features/name")
    h5write(rep("Gene Expression", length(features)), h5file, "matrix/features/feature_type")
    h5write(rep("GRCh38", length(features)), h5file, "matrix/features/genome")
}

prefix <- "/lustre/scratch127/gengen/projects_v2/dual_tf/MORF_10/Multiome/cellranger-arc/cellranger-arc202_count_"
atacFiles <- c("morf_b1" = paste0(prefix, "e114917b4fbd8ec64708b8626edcca63/atac_fragments.tsv.gz"),
               "morf_b2" = paste0(prefix, "17c16ebd84bcf0d3a5b9c32235070443/atac_fragments.tsv.gz"),
               "morf_b3" = paste0(prefix, "d3bfc0670591f217496be4a028c6bfeb/atac_fragments.tsv.gz"),
               "morf_b4" = paste0(prefix, "78ab1a8f128d665f687011190ad63365/atac_fragments.tsv.gz"),
               "morf_b5" = paste0(prefix, "d3182970f44fc9fedbf4d84385c9470f/atac_fragments.tsv.gz"),
               "morf_b6" = paste0(prefix, "4dea0ff7eeae5dd3d8c156c404de6f9f/atac_fragments.tsv.gz"),
               "morf_b7" = paste0(prefix, "8d7cec729ebd63f6bffeb3aecf6ef800/atac_fragments.tsv.gz"),
               "morf_b8" = paste0(prefix, "6027535ef50f68d62e0ce589da9cb022/atac_fragments.tsv.gz"))
# rnaFiles <- c("morf_b1" = paste0(prefix, "e114917b4fbd8ec64708b8626edcca63/filtered_feature_bc_matrix.h5"),
#               "morf_b2" = paste0(prefix, "17c16ebd84bcf0d3a5b9c32235070443/filtered_feature_bc_matrix.h5"),
#               "morf_b3" = paste0(prefix, "d3bfc0670591f217496be4a028c6bfeb/filtered_feature_bc_matrix.h5"),
#               "morf_b4" = paste0(prefix, "78ab1a8f128d665f687011190ad63365/filtered_feature_bc_matrix.h5"),
#               "morf_b5" = paste0(prefix, "d3182970f44fc9fedbf4d84385c9470f/filtered_feature_bc_matrix.h5"),
#               "morf_b6" = paste0(prefix, "4dea0ff7eeae5dd3d8c156c404de6f9f/filtered_feature_bc_matrix.h5"),
#               "morf_b7" = paste0(prefix, "8d7cec729ebd63f6bffeb3aecf6ef800/filtered_feature_bc_matrix.h5"),
#               "morf_b8" = paste0(prefix, "6027535ef50f68d62e0ce589da9cb022/filtered_feature_bc_matrix.h5"))
prefix <- "/lustre/scratch127/gengen/projects_v2/dual_tf/MORF_10/Multiome/analysis_by_hgi/4_archr/qced_h5"
rnaFiles <- c("morf_b1" = paste0(prefix, "/morf10_b1_qced.h5"),
              "morf_b2" = paste0(prefix, "/morf10_b2_qced.h5"),
              "morf_b3" = paste0(prefix, "/morf10_b3_qced.h5"),
              "morf_b4" = paste0(prefix, "/morf10_b4_qced.h5"),
              "morf_b5" = paste0(prefix, "/morf10_b5_qced.h5"),
              "morf_b6" = paste0(prefix, "/morf10_b6_qced.h5"),
              "morf_b7" = paste0(prefix, "/morf10_b7_qced.h5"),
              "morf_b8" = paste0(prefix, "/morf10_b8_qced.h5"))


ArrowFiles <- createArrowFiles(inputFiles = atacFiles,
                               sampleNames = names(atacFiles),
                               minTSS = 4,
                               minFrags = 1000,
                               TileMatParams = list(tileSize = 1000),
                               addTileMat = TRUE,
                               addGeneScoreMat = TRUE,
                               threads = 16)
proj_morf10_atac <- ArchRProject(ArrowFiles = ArrowFiles)

seRNA <- import10xFeatureMatrix(input = rnaFiles,
                                names = names(rnaFiles),
                                strictMatch = TRUE)

