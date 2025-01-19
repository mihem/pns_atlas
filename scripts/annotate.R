#===============================================================================
# Cell Type Annotation Script
#===============================================================================
# Purpose: Annotate cell clusters and generate visualization plots, including:
# - Loading and applying cluster annotations from external files
# - Setting up custom visualization parameters
# - Generating UMAP plots with different groupings
# - Creating feature plots for marker genes
# - Generating dot plots for cell type markers
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
library(scMisc)
library(Polychrome)
library(readxl)

# general settings  ----
options(warn = 0)
options(Seurat.object.assay.version = "v5")
future::plan("multicore", workers = 6)
options(future.globals.maxSize = 16000 * 1024^2)
conflicts_prefer(base::setdiff)
conflicts_prefer(base::as.data.frame)
my_cols_25 <- pals::cols25()
my_cols_50 <- unname(Polychrome::createPalette(50, pals::cols25()))
markers_pns <- read_csv("/home/mischko/Documents/beruf/forschung/markers/markers_pns.csv")

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)

# annotate clusters ----
# Apply cluster annotations from external file and set up visualization parameters
Idents(sc_merge) <- sc_merge$RNA_snn_res.0.7

annotations <- readxl::read_xlsx(file.path("lookup", "cluster_annotation.xlsx")) |>
    mutate(cluster = as.character(cluster)) |>
    arrange(cluster)

names(annotations$final) <- levels(sc_merge)
sc_merge <- RenameIdents(sc_merge, annotations$final)
 
# save custom cluster and condition order and colors in seurat object
cluster_order <- readxl::read_xlsx(file.path("lookup", "cluster_order.xlsx"))
sc_merge@misc$cluster_order <- cluster_order[!is.na(cluster_order)]
sc_merge@misc$cluster_col <- setNames(pals::cols25(length(sc_merge@misc$cluster_order)), sc_merge@misc$cluster_order)
sc_merge@misc$level2_order <- c("CTRL", "VN", "CIDP", "CIAP", "PPN", "DPN", "OIN", "ONIN")
sc_merge@misc$level2_cols <- setNames(pals::cols25(length(sc_merge@misc$level2_order)), sc_merge@misc$level2_order)

#sanity check
setdiff(sc_merge@misc$cluster_order, levels(sc_merge))
setdiff(levels(sc_merge), sc_merge@misc$cluster_order)
any(duplicated(sc_merge@misc$cluster_order))

# save ordered annotations in meta.data
sc_merge$cluster <- factor(Idents(sc_merge), levels = sc_merge@misc$cluster_order)
Idents(sc_merge) <- sc_merge$cluster

DefaultAssay(sc_merge) <- "RNA"

# Generate main UMAP visualization with cluster annotations
umap <-
    DimPlot(sc_merge, reduction = "umap.scvi.full", pt.size = .1, raster = FALSE, alpha = 0.1, group.by = "cluster", cols = sc_merge@misc$cluster_col, label = FALSE) +
    theme_rect() +
    NoLegend() + 
    xlab("") + 
    ylab("") + 
    ggtitle("")

ggsave(plot = umap, file.path("results", "umap", "scvi_umap_full_annotated.png"), width = 10, height = 8)
ggsave(plot = umap, file.path("results", "umap", "scvi_umap_full_annotated_no_axis.png"), width = 10, height = 8)
ggsave(plot = umap, file.path("results", "umap", "scvi_umap_full_annotated.pdf"), width = 10, height = 8)

# final UMAP plot by center
umap_center <-
    DimPlot(sc_merge, reduction = "umap.scvi.full", pt.size = .1, raster = FALSE, alpha = 0.1, group.by = "center", cols = pals::cols25(3), label = FALSE) +
    theme_rect() +
    xlab("UMAP1") + 
    ylab("UMAP2")

ggsave(plot = umap_center, file.path("results", "umap", "scvi_umap_group_center_full.png"), width = 10, height = 8)
ggsave(plot = umap_center, file.path("results", "umap", "scvi_umap_group_center_full.pdf"), width = 10, height = 8)


# final UMAP plot by sample
umap_sample <-
    DimPlot(sc_merge, reduction = "umap.scvi.full", pt.size = .1, raster = FALSE, alpha = 0.1, group.by = "sample", cols = my_cols_50, label = FALSE) +
    theme_rect() +
    xlab("UMAP1") + 
    ylab("UMAP2")

ggsave(plot = umap_sample, file.path("results", "umap", "scvi_umap_group_sample_full.png"), width = 10, height = 8)
ggsave(plot = umap_sample, file.path("results", "umap", "scvi_umap_group_sample_full.pdf"), width = 10, height = 8)

# save object ---
qs::qsave(sc_merge, file.path("objects", "sc_merge.qs"))

# feature plots ----
lapply(unique(markers_pns$cell_source),
  FUN = purrr::possibly(scMisc::fPlotCustom),
  object = sc_merge,
  markers = markers_pns,
  reduction = "umap.scvi.full"
)

# dot plots ----
DefaultAssay(sc_merge) <- "RNA"

dotPlot(
  path = file.path("lookup", "markers.csv"),
  object = sc_merge,
  par = "dotplot_long",
  dot_min = 0.01,
  height = 8,
  width = 26
)

dotPlot(
  path = file.path("lookup", "markers.csv"),
  object = sc_merge,
  par = "dotplot_short",
  dot_min = 0.01,
  height = 7,
  width = 16
)

dotPlot(
  path = file.path("lookup", "markers.csv"),
  object = sc_merge,
  par = "Macro18",
  dot_min = 0.01,
  height = 7.5,
  width = 20
)

Idents(sc_merge) <- factor(sc_merge$cluster, levels = rev(sc_merge@misc$cluster_order))

dotPlot(
  path = file.path("lookup", "markers.csv"),
  object = sc_merge,
  par = "dotplot_jolien",
  dot_min = 0.01,
  height = 7,
  width = 15
)

# specialized dot plots for specific cell types ----
ec <- subset(sc_merge, subset = RNA_snn_res.0.7 %in% c("7", "10", "11", "19"))

lapply(
  c("ec_bnb_jolien", "EC_BNB", "Yang_2022_venous", "Yang_2022_capillary", "Yang_2022_arterial"),
  FUN = function(x) {
    dotPlot(
      path = file.path("lookup", "markers.csv"),
      object = ec,
      par = x,
      dot_min = 0.01,
      height = 4,
      width = 10
    )
  }
)

lapply(
  c("Garcia_2022_capillary", "Garcia_2022_arterial", "Garcia_2022_venous"),
  FUN = function(x) {
    dotPlot(
      path = file.path("lookup", "markers.csv"),
      object = ec,
      par = x,
      dot_min = 0.01,
      height = 4,
      width = 10
    )
  }
)

vsmc_pc <- subset(sc_merge, subset = RNA_snn_res.0.7 %in% c("3", "6", "17"))

lapply(
  c("telocytes_jolien"),
  FUN = function(x) {
    dotPlot(
      path = file.path("lookup", "markers.csv"),
      object = vsmc_pc,
      par = x,
      dot_min = 0.01,
      height = 4,
      width = 10
    )
  }
)
