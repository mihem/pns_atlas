#===============================================================================
# Data Integration Script
#===============================================================================
# Purpose: Integrate single-cell data across samples and map to reference datasets
#
# Key Analysis Steps:
# 1. Data integration using scVI
# 2. Initial analysis on sketched subset for efficiency
# 3. Project integration to full dataset
# 4. Reference mapping to multiple datasets:
#    - PNS atlas (Milbrandt)
#    - EC atlas
#    - Suter sciatic nerve atlas (P1/P60)
#    - BBB/vascular references
#    - ROSMAP dataset
#===============================================================================

# Load libraries and set up environment ----
library(Seurat)
library(BPCells)
library(SeuratObject)
library(SeuratWrappers)
library(tidyverse)
library(writexl)
library(patchwork)
library(conflicted)
library(qs)
library(pals)
library(homologene)
library(scMisc)
library(pals)
library(Polychrome)

# general settings  ----
options(warn = 0)
options(Seurat.object.assay.version = "v5")
future::plan("multicore", workers = 6)
options(future.globals.maxSize = 16000 * 1024^2)
conflicts_prefer(base::setdiff)
my_cols_25 <- pals::cols25()
my_cols_50 <- unname(Polychrome::createPalette(50, pals::cols25()))

# Load and normalize data ----
sc_merge_pre <- qs::qread(file.path("objects", "sc_merge_pre.qs"), nthread = 4)

# normalize ----
sc_merge <- NormalizeData(sc_merge_pre, verbose = TRUE, normalization.method = "LogNormalize", scale.factor = 10000)
sc_merge <- FindVariableFeatures(sc_merge, selection.method = "vst", nfeatures = 2000)

# Sketch data for initial analysis ----
# Create efficient subset using LeverageScore sampling
sc_merge <- SketchData(object = sc_merge, ncells = 5000, method = "LeverageScore", sketched.assay = "sketch")

## DefaultAssay(sc_merge) <- "RNA"
DefaultAssay(sc_merge) <- "sketch"

sc_merge <-
  sc_merge |>
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) |>
  ScaleData() |>
  RunPCA()

# Integrate data using scVI ----
# Note: Requires scvi-tools installation in conda env
# comment: requires scvi installation, problematic with radian (which forces a specific conda env), and with scvi-tools 1.0.2
set.seed(123)
sc_merge <- IntegrateLayers(
  object = sc_merge,
  method = scVIIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.scvi",
  conda_env = "~/miniconda3/envs/scvi",
  group.by = "sample",
  verbose = TRUE
)

# run umap on sketch integration
sc_merge <- RunUMAP(sc_merge, reduction = "integrated.scvi", reduction.name = "umap.scvi.sketch", dims = 1:30)

# plot integration
umap_group_sample_sketch <- DimPlot(sc_merge, reduction = "umap.scvi.sketch", group.by = "sample", raster = FALSE, pt.size = .1, alpha = .1, cols = my_cols_50) + theme_rect()
ggsave(plot = umap_group_sample_sketch, file.path("results", "umap", "scvi_umap_group_sample_sketch.png"), width = 10, height = 8)

umap_group_center_sketch <- DimPlot(sc_merge, reduction = "umap.scvi.sketch", group.by = "center", raster = FALSE, pt.size = .1, alpha = .1, cols = my_cols_25) + theme_rect()
ggsave(plot = umap_group_center_sketch, file.path("results", "umap", "scvi_umap_group_center_sketch.png"), width = 10, height = 8)

qs::qsave(sc_merge, file.path("objects", "sc_merge.qs"))

# Map to reference datasets ----
# Helper functions for mapping
convertRownames <- function(seu_object) {
  lookup <- homologene::mouse2human(rownames(seu_object), db = homologene::homologeneData2)
  new_rownames <- lookup$humanGene[match(rownames(seu_object), lookup$mouseGene)]
  rownames(seu_object@assays$RNA@counts) <- new_rownames
  rownames(seu_object@assays$RNA@data) <- new_rownames
  # rownames(seu_object@assays$RNA@scale.data) <- new_rownames
  #remove columns with NA
  features_keep <- rownames(seu_object)[!is.na(rownames(seu_object))]
  obj_new <- subset(seu_object, features = features_keep)
  rownames(obj_new@assays$RNA@meta.features) <- rownames(obj_new)
  return(obj_new)
}

mapSeurat <- function(ref, query) {
  reference_list <- list(ref = ref, query = query)
  features <- SelectIntegrationFeatures(object.list = reference_list)
  anchors <- FindTransferAnchors(
    reference = reference_list$ref,
    query = reference_list$query,
    normalization.method = "LogNormalize",
    features = features
  )
  predictions <- TransferData(anchorset = anchors, refdata = reference_list$ref$cluster)
  return(predictions)
}

storePred <- function(predictions, label_col, score_col, seu_obj) {
  predictions_prep <-
    predictions |>
    tibble::rownames_to_column("barcode") |>
    dplyr::select(predicted.id, prediction.score.max, barcode) |>
    dplyr::mutate(predicted.id = ifelse(prediction.score.max < 0.3, "unknown", predicted.id)) |>
    tibble::as_tibble() |>
    dplyr::rename(
      {{ label_col }} := predicted.id,
      {{ score_col }} := prediction.score.max
    )

  seu_obj@meta.data <-
    seu_obj@meta.data |>
    tibble::rownames_to_column("barcode") |>
    dplyr::left_join(predictions_prep, by = "barcode") |>
    tibble::column_to_rownames(var = "barcode")

  return(seu_obj)
}

# Load and prepare reference datasets ----
pns_sn_sciatic_milbrandt <- qs::qread("/home/mischko/Documents/beruf/forschung/scRNA_reference/pns_atlas_milbrandt/pns_sn_sciatic_GSE182098.qs", nthreads = 4)
DimPlot(pns_sn_sciatic_milbrandt, label = TRUE)
dplyr::count(pns_sn_sciatic_milbrandt@meta.data, cluster)

suter_p1 <- readRDS("/home/mischko/Documents/beruf/forschung/scRNA_reference/sciatic_nerve_atlas_suter/10xGeno_P1.rds")
suter_p1$cluster <- Idents(suter_p1)

suter_p60 <- readRDS("/home/mischko/Documents/beruf/forschung/scRNA_reference/sciatic_nerve_atlas_suter/10xGeno_P60.rds")

suter_p1_diet <- DietSeurat(suter_p1)
suter_p60_diet <- DietSeurat(suter_p60)

suter_merge <- merge(x = suter_p1_diet, y = suter_p60_diet)
dplyr::count(suter_merge@meta.data, cluster)

ec_atlas <- qs::qread("/home/mischko/Documents/beruf/forschung/scRNA_reference/ec_atlas/ec_atlas.qs")
ts_ec_small <- qs::qread("/home/mischko/Documents/beruf/forschung/scRNA_reference/tabula_sapiens/ts_ec_small.qs")
bbb_vascular <- qs::qread("/home/mischko/Documents/beruf/forschung/scRNA_reference/bbb/bbb_vascular.qs")
rosmap <- qs::qread("/home/mischko/Documents/beruf/forschung/scRNA_reference/bbb/rosmap.qs")

# prepare reference
human_pns_sciatic_milbrandt <- convertRownames(pns_sn_sciatic_milbrandt)
human_ec_atlas <- convertRownames(ec_atlas)
human_suter_p60 <- convertRownames(suter_p60)
human_suter_p1 <- convertRownames(suter_p1)
human_suter_merge <- convertRownames(suter_merge)

# Perform reference mapping ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"))
sc_merge[["sketch"]] <- JoinLayers(sc_merge[["sketch"]])

DefaultAssay(sc_merge) <- "sketch"

predictions_milbrandt <- mapSeurat(ref = human_pns_sciatic_milbrandt, query = sc_merge)
predictions_ec <- mapSeurat(ref = human_ec_atlas, query = sc_merge)
predictions_ts_ec <- mapSeurat(ref = ts_ec_small, query = sc_merge)
predictions_suter_p1 <- mapSeurat(ref = human_suter_p1, query = sc_merge)
predictions_suter_p60 <- mapSeurat(ref = human_suter_p60, query = sc_merge)
predictions_suter_merge <- mapSeurat(ref = human_suter_merge, query = sc_merge)
predictions_bbb <- mapSeurat(ref = bbb_vascular, query = sc_merge)
predictions_rosmap <- mapSeurat(ref = rosmap, query = sc_merge)


sc_merge <- storePred(predictions_milbrandt, label_col = "milbrandt_sciatic_label", score_col = "milbrandt_sciatic_score", seu_obj = sc_merge)
sc_merge <- storePred(predictions_ec, label_col = "ec_atlas_label", score_col = "ec_atlas_score", seu_obj = sc_merge)
sc_merge <- storePred(predictions_ts_ec, label_col = "ts_ec_label", score_col = "ts_ec_score", seu_obj = sc_merge)
sc_merge <- storePred(predictions_suter_p1, label_col = "suter_p1_label", score_col = "suter_p1_score", seu_obj = sc_merge)
sc_merge <- storePred(predictions_suter_p60, label_col = "suter_p60_label", score_col = "suter_p60_score", seu_obj = sc_merge)
sc_merge <- storePred(predictions_suter_merge, label_col = "suter_merge_label", score_col = "suter_merge_score", seu_obj = sc_merge)
sc_merge <- storePred(predictions_bbb, label_col = "bbb_label", score_col = "bbb_score", seu_obj = sc_merge)
sc_merge <- storePred(predictions_rosmap, label_col = "rosmap_label", score_col = "rosmap_score", seu_obj = sc_merge)

#sanity check
table(sc_merge@meta.data$heming_sural_label)
table(sc_merge@meta.data$milbrandt_sciatic_label)
table(sc_merge@meta.data$ec_atlas_label)
table(sc_merge@meta.data$ts_ec_label)
table(sc_merge@meta.data$suter_p1_label)
table(sc_merge@meta.data$suter_p60_label)
table(sc_merge@meta.data$bbb_label)

pred_plot_heming_sketch <-
  DimPlot(sc_merge, reduction = "umap.scvi.sketch", group.by = "heming_sural_label", raster = FALSE, pt.size = .1, alpha = .1, cols = my_cols_25, label = TRUE) +
  theme_rect() +
  NoLegend()
ggsave(plot = pred_plot_heming_sketch, file.path("results", "map", "map_heming_sketch.png"), width = 8, height = 8)

pred_plot_milbrandt_sketch <-
 DimPlot(sc_merge, reduction = "umap.scvi.sketch", group.by = "milbrandt_sciatic_label", raster = FALSE, pt.size = .1, alpha = .1, cols = my_cols_25, label = TRUE) +
  theme_rect() + 
  NoLegend()
ggsave(plot = pred_plot_milbrandt_sketch, file.path("results", "map", "map_milbrandt_sketch.png"), width = 8, height = 8)

pred_plot_ec_atlas_sketch <-
 DimPlot(sc_merge, reduction = "umap.scvi.sketch", group.by = "ec_atlas_label", raster = FALSE, pt.size = .1, alpha = .1, cols = my_cols_25, label = TRUE) +
  theme_rect() + 
  NoLegend()
ggsave(plot = pred_plot_ec_atlas_sketch, file.path("results", "map", "map_ec_atlas_sketch.png"), width = 8, height = 8)

pred_plot_ec_atlas_sketch <-
 DimPlot(sc_merge, reduction = "umap.scvi.sketch", group.by = "suter_p60_label", raster = FALSE, pt.size = .1, alpha = .1, cols = my_cols_25, label = TRUE) +
  theme_rect() + 
  NoLegend()

ggsave(plot = pred_plot_ec_atlas_sketch, file.path("results", "map", "map_ec_atlas_sketch.png"), width = 8, height = 8)

# Project sketched integration onto full dataset ---
DefaultAssay(sc_merge) <- "sketch"
sc_merge[["sketch"]] <- split(sc_merge[["sketch"]], f = sc_merge$sample)

sc_merge <- ProjectIntegration(
  object = sc_merge,
  sketched.assay = "sketch",
  assay = "RNA",
  reduction = "integrated.scvi")

sc_merge <- ProjectData(
  object = sc_merge,
  assay = "RNA",
  sketched.assay = "sketch",
  sketched.reduction = "integrated.scvi",
  full.reduction = "integrated.scvi.full",
  dims = 1:30,
  refdata = list(
    heming_sural_label_full = "heming_sural_label",
    milbrandt_sciatic_label_full = "milbrandt_sciatic_label",
    ec_atlas_label_full = "ec_atlas_label",
    ts_ec_label_full = "ts_ec_label",
    suter_p1_label_full = "suter_p1_label",
    suter_p60_label_full = "suter_p60_label",
    suter_merge_label_full = "suter_merge_label"
  )
)

# run umap based on full integrated scvi
sc_merge <- RunUMAP(
  object = sc_merge,
  reduction = "integrated.scvi.full",
  reduction.name = "umap.scvi.full",
  dims = 1:30)

# now that we have projected the full dataset, switch back to analyzing all cells
DefaultAssay(sc_merge) <- "RNA"
qs::qsave(sc_merge, file.path("objects", "sc_merge.qs"))

# Generate visualizations ----
# Plot UMAPs and prediction results
umap_group_sample <-
  DimPlot(sc_merge, reduction = "umap.scvi.full", pt.size = .1, raster = FALSE, alpha = 0.1, group.by = "sample", cols = my_cols_50) +
  theme_rect() 
ggsave(plot = umap_group_sample, file.path("results", "umap", "scvi_umap_group_sample_full.png"), width = 10, height = 8)

umap_group_center <-
  DimPlot(sc_merge, reduction = "umap.scvi.full", pt.size = .1, raster = FALSE, alpha = 0.1, group.by = "center", cols = my_cols_50) +
  theme_rect()
ggsave(plot = umap_group_center, file.path("results", "umap", "scvi_umap_group_center_full.png"), width = 10, height = 8)

pred_plot_heming_full <-
  DimPlot(sc_merge, reduction = "umap.scvi.full", group.by = "heming_sural_label_full", raster = FALSE, pt.size = .1, alpha = .1, cols = my_cols_25, label = TRUE) +
  theme_rect() +
  NoLegend()
ggsave(plot = pred_plot_heming_full, file.path("results", "map", "map_heming_full.png"), width = 8, height = 8)

pred_plot_milbrandt_full <-
  DimPlot(sc_merge, reduction = "umap.scvi.full", group.by = "milbrandt_sciatic_label_full", raster = FALSE, pt.size = .1, alpha = .1, cols = my_cols_25, label = FALSE, label.size = 6) +
  theme_rect() +
  NoLegend() +
  xlab("UMAP1") +
  ylab("UMAP2") +
  ggtitle("Yim et al.")

ggsave(plot = pred_plot_milbrandt_full, file.path("results", "map", "map_milbrandt_full.png"), width = 8, height = 8)

pred_plot_suter_p1_full <-
 DimPlot(sc_merge, reduction = "umap.scvi.full", group.by = "suter_p1_label_full", raster = FALSE, pt.size = .1, alpha = .1, cols = my_cols_25, label = TRUE) +
  theme_rect() + 
  NoLegend() + 
ggsave(plot = pred_plot_suter_p1_full, file.path("results", "map", "map_suter_p1_full.png"), width = 8, height = 8)

pred_plot_suter_p60_full <-
 DimPlot(sc_merge, reduction = "umap.scvi.full", group.by = "suter_p60_label_full", raster = FALSE, pt.size = .1, alpha = .1, cols = my_cols_25, label = FALSE, label.size = 6) +
  theme_rect() + 
  NoLegend() + 
  xlab("UMAP1") +
  ylab("UMAP2") + 
  ggtitle("Gerber et al. p60")

ggsave(plot = pred_plot_suter_p60_full, file.path("results", "map", "map_suter_p60_full.png"), width = 8, height = 8)

pred_plot_suter_merge_full <-
 DimPlot(sc_merge, reduction = "umap.scvi.full", group.by = "suter_merge_label_full", raster = FALSE, pt.size = .1, alpha = .1, cols = my_cols_25, label = TRUE) +
  theme_rect() 
ggsave(plot = pred_plot_suter_merge_full, file.path("results", "map", "map_suter_merge_full.png"), width = 8, height = 8)

# this can only be run after cluster.R was run
ec <- subset(sc_merge, subset = RNA_snn_res.0.7 %in% c("7", "10", "11", "19"))

pred_plot_ec_atlas_full <-
 DimPlot(ec, reduction = "umap.scvi.full", group.by = "ec_atlas_label_full", raster = FALSE, pt.size = .1, alpha = .1, cols = my_cols_25, label = FALSE) +
  theme_rect()
ggsave(plot = pred_plot_ec_atlas_full, file.path("results", "map", "map_ec_atlas_full.png"), width = 8, height = 8)

pred_plot_ts_ec <-
 DimPlot(ec, reduction = "umap.scvi.full", group.by = "ts_ec_label_full", raster = FALSE, pt.size = .1, alpha = .1, cols = my_cols_25, label = FALSE) +
  theme_rect() 
ggsave(plot = pred_plot_ts_ec, file.path("results", "map", "map_ts_ec.png"), width = 8, height = 8)

sc <- subset(sc_merge, subset = cluster %in% c("mySC", "nmSC", "repairSC", "damageSC"))

pred_plot_sc_suter_p1_full <-
 DimPlot(sc, reduction = "umap.scvi.full", group.by = "suter_p1_label_full", raster = FALSE, pt.size = .1, alpha = .1, cols = my_cols_25, label = FALSE) +
  theme_rect()
ggsave(plot = pred_plot_sc_suter_p1_full, file.path("results", "map", "map_sc_suter_p1_full.png"), width = 8, height = 8)

pred_plot_bbb <-
 DimPlot(sc_merge, reduction = "umap.scvi.full", group.by = "bbb_label", raster = FALSE, pt.size = .1, alpha = .1, cols = my_cols_25, label = FALSE) +
  theme_rect()
ggsave(plot = pred_plot_bbb, file.path("results", "map", "map_bbb_full.png"), width = 8, height = 8)

# only keep large clusters
ec_rosmap <- subset(ec, subset = rosmap_label %in% c("aEndo", "capEndo", "vEndo"))
table(ec_rosmap$rosmap_label)

pred_plot_rosmap <-
 DimPlot(ec_rosmap, reduction = "umap.scvi.full", group.by = "rosmap_label", raster = FALSE, pt.size = .1, alpha = .1, cols = my_cols_25, label = FALSE)  + 
  theme_rect() +
  NoLegend() + 
  xlim(-12, -5) +
  ylim(2, 10) +
  xlab("UMAP1") +
  ylab("UMAP2") + 
  ggtitle("ROSMAP vascular cells")

ggsave(plot = pred_plot_rosmap, file.path("results", "map", "map_rosmap_full.png"), width = 4, height = 4)

