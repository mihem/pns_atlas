# libraries
library(Seurat)
library(scSHC)
library(tidyverse)
library(scMisc)

# general settings  ----
options(warn = 0)
options(Seurat.object.assay.version = "v5")
my_cols_50 <- unname(Polychrome::createPalette(50, pals::cols25()))

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)

# downsample subset for shc
Idents(sc_merge) <- sc_merge$RNA_snn_res.0.7

sc_small <- subset(sc_merge, downsample = 1000)

qs::qsave(sc_small, file.path("objects", "sc_small.qs"))

counts_small <- as(object = sc_small[["RNA"]]$counts, Class = "dgCMatrix")
qs::qsave(counts_small, file.path("objects", "counts_small.qs"))
qs::qsave(sc_small@meta.data, file.path("objects", "sc_small_meta_data.qs"))

# reload because of high memory usage
counts_small <- qs::qread(file.path("objects", "counts_small.qs"))
meta_data <- qs::qread(file.path("objects", "sc_small_meta_data.qs"))

new_clusters <-
    scSHC::testClusters(
        counts_small,
        cluster_ids = as.character(meta_data$RNA_snn_res.0.7),
        batch = meta_data$sample,
        parallel = FALSE,
        cores = 1
    )

qs::qsave(new_clusters, file.path("objects", "sc_shc_clusters.qs"))

new_clusters <- qs::qread(file.path("objects", "sc_shc_clusters.qs"))
sc_small <- qs::qread(file.path("objects", "sc_small.qs"))

str(new_clusters[[1]])
sc_small

# Assign SHC clusters to Seurat object
sc_small$shc <- new_clusters[[1]]
sc_small$shc <- gsub(x = sc_small$shc, pattern = "new", replacement = "cl")

# Generate a UMAP plot with shc clusters
umap_shc <- DimPlot(sc_small, reduction = "umap.scvi.full", pt.size = .5, raster = FALSE, alpha = 0.1, group.by = "shc", cols = my_cols_50, label = TRUE) +
    theme_rect() +
    NoLegend()

ggsave(plot = umap_shc, file.path("results", "umap", "scvi_umap_small_shc.png"), width = 10, height = 8)

qs::qsave(sc_small, file.path("objects", "sc_small.qs"))

dplyr::count(sc_small@meta.data, shc)

new_clusters
