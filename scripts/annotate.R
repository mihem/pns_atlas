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

remotes::install_github("mihem/scMisc")
detach(package:scMisc, unload = TRUE)

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

# final annotated UMAP plot by cluster
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

# dotplot ---
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
  par = "dotplot_jolien",
  dot_min = 0.01,
  height = 7,
  width = 15
)