#===============================================================================
# Clustering Analysis Script
#===============================================================================
# Purpose: Process single-cell RNA sequencing data for cluster analysis:
# - Perform clustering at multiple resolutions (0.2-0.8)
# - Generate UMAP visualizations of clustering results
# - Identify cluster-specific marker genes
# - Export results for downstream analysis
#===============================================================================

# libraries  ----
library(Seurat)
library(BPCells)
library(SeuratObject)
library(tidyverse)
library(writexl)
library(patchwork)
library(conflicted)
library(qs)
library(pals)
library(Polychrome)
library(scMisc)

# general settings  ----
options(warn = 0)
options(Seurat.object.assay.version = "v5")
future::plan("multicore", workers = 6)
options(future.globals.maxSize = 16000 * 1024^2)
conflicts_prefer(base::setdiff)
my_cols_25 <- pals::cols25()
my_cols_50 <- unname(Polychrome::createPalette(50, pals::cols25()))
markers_pns <- read_csv("/home/mischko/Documents/beruf/forschung/markers/markers_pns.csv")

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)

str(sc_merge@meta.data)

# find clusters ----
DefaultAssay(sc_merge) <- "RNA"

sc_merge <- FindNeighbors(sc_merge, reduction = "integrated.scvi.full", dims = 1:30, assay = "RNA")

for (res in seq(from = 0.2, to = 0.8, by = 0.1)) {
  sc_merge <- FindClusters(sc_merge, resolution = res)
}

qs::qsave(sc_merge, file.path("objects", "sc_merge.qs"))

# plot clustering ----
umap <-
    DimPlot(sc_merge, reduction = "umap.scvi.full", pt.size = .1, raster = FALSE, alpha = 0.1, group.by = "RNA_snn_res.0.7", cols = my_cols_50, label = TRUE) +
    theme_rect() +
    NoLegend()
ggsave(plot = umap, file.path("results", "umap", "scvi_umap_full.png"), width = 10, height = 8)

resolutions <- paste0("RNA_snn_res.", seq(from = 0.2, to = 0.8, by = 0.1))

umap_list <- list()
for (res in resolutions) {
    umap_list[[res]] <-
        DimPlot(sc_merge, reduction = "umap.scvi.full", pt.size = .1, raster = FALSE, alpha = 0.1, group.by = res, cols = my_cols_50, label = TRUE) +
        theme_rect() +
        NoLegend()
}

umap_list <- patchwork::wrap_plots(umap_list, ncol = 2)
ggsave(plot = umap_list, file.path("results", "umap", "scvi_umap_full_resolutions.png"), width = 16, height = 30)

# find markers ----
sc_merge <- JoinLayers(sc_merge)

# find markers helper function
findMarkers <- function(ident1, ident2 = NULL, object, only_pos, min_pct, logfc_threshold, assay = assay) {
  result <- Seurat::FindMarkers(object, ident.1 = ident1, ident.2 = ident2, min.pct = min_pct, logfc.threshold = logfc_threshold, only.pos = only_pos, assay = assay) |>
    tibble::rownames_to_column("gene") |>
    dplyr::filter(p_val_adj < 0.05) |>
    dplyr::relocate(gene, avg_log2FC, p_val, p_val_adj) |>
    dplyr::arrange(desc(avg_log2FC))
  return(result)
}

topmarkers <-
  lapply(
    unique(sc_merge@misc$cluster_order),
    function(x) {
      message("Processing cluster ", x)
      try(findMarkers(ident1 = x, object = sc_merge, only_pos = TRUE, min_pct = 0.1, logfc_threshold = 0.25, assay = "RNA"))
    }
  )

names(topmarkers) <- sc_merge@misc$cluster_order

writexl::write_xlsx(topmarkers, file.path("results", "de", "topmarkers_final.xlsx"))
