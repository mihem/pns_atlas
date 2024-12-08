# ===============================================================================
# Immune Cell Analysis Script
# ===============================================================================
# Purpose: Characterize and annotate immune cell populations from snRNA-seq data
#
# Key Analysis Steps:
# 1. Subset immune cells from main dataset
# 2. Reference mapping using multiple references:
# 3. Integration and dimensional reduction
# 4. Clustering optimization
# 5. Manual annotation based on markers and reference mapping
#
# ===============================================================================

# Load required libraries and set up environment ----
library(Seurat)
library(BPCells)
library(SeuratObject)
library(tidyverse)
library(patchwork)
library(conflicted)
library(qs)
library(scMisc)
library(Azimuth)
library(SeuratData)
library(SeuratWrappers)

# general settings  ----
options(warn = 0)
options(Seurat.object.assay.version = "v5")
future::plan("multicore", workers = 6)
options(future.globals.maxSize = 16000 * 1024^2)
conflicts_prefer(base::setdiff)
my_cols_25 <- pals::cols25()
my_cols_50 <- unname(Polychrome::createPalette(50, pals::cols25()))

# Load and subset immune cell populations ----
# load preprocessed data ----
ic <- qread(file.path("objects", "ic.qs"), nthread = 4)

# subset immune cells and recluster ----
# Extract immune cell clusters and recluster
ic <- subset(sc_merge, cluster %in% c("Macro1", "Macro2", "Granulo", "B", "T_NK", "Mast"))
ic$cluster <- droplevels(ic$cluster)

# Reference mapping with Azimuth PBMC reference ----
# run Azimuth mapping to references ----
# Map to PBMC reference
ic_azimuth <- RunAzimuth(query = ic, reference = "pbmcref", assay = "RNA")
qsave(ic_azimuth, file.path("objects", "ic_azimuth.qs"))

# prepare predictions ----
predictions_pbmcref <-
  data.frame(
    ic_level1 = ic_azimuth$predicted.celltype.l1,
    ic_level1_score = ic_azimuth$predicted.celltype.l1.score,
    ic_level2 = ic_azimuth$predicted.celltype.l2,
    ic_level2_score = ic_azimuth$predicted.celltype.l2.score,
    ic_level3 = ic_azimuth$predicted.celltype.l3,
    ic_level3_score = ic_azimuth$predicted.celltype.l3.score
  ) |>
  mutate(ic_level1 = ifelse(ic_level1_score < 0.4, "unknown", ic_level1)) |>
  mutate(ic_level2 = ifelse(ic_level2_score < 0.4, "unknown", ic_level2)) |>
  mutate(ic_level3 = ifelse(ic_level3_score < 0.4, "unknown", ic_level3))

# add to seurat object ---
ic <- AddMetaData(ic, predictions_pbmcref)

# plot predictions ----
pred_plot_ic_lev2 <-
  DimPlot(ic, reduction = "umap.scvi.full", group.by = "ic_level2", raster = FALSE, pt.size = .1, alpha = .1, cols = rev(my_cols_25), label = FALSE) +
  theme_rect()

ggsave(plot = pred_plot_ic_lev2, file.path("results", "map", "map_ic_lev2.png"), width = 15, height = 10)

# Optimize data for integration ----
ic_slim <- DietSeurat(
  ic,
  assay = c("RNA"),
  layers = c("counts"),
)

# convert to sparse matrices ----
ic_slim$RNA$counts <- as(object = ic_slim$RNA$counts, Class = "dgCMatrix")
ic_slim$RNA$data <- as(object = ic_slim$RNA$data, Class = "dgCMatrix")
ic_slim[["RNA"]]$scale.data <- NULL

# Perform multi-modal integration ----
ic_slim[["RNA"]] <- split(
  x = ic_slim[["RNA"]],
  f = ic_slim$sample
)

ic_slim <-
  ic_slim |>
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) |>
  ScaleData() |>
  RunPCA()

ic_slim <- IntegrateLayers(
  object = ic_slim,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = TRUE,
)

ic_slim <- IntegrateLayers(
  object = ic_slim,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "rpca",
  verbose = TRUE,
  k.weight = 30
)

set.seed(123)
ic_slim <- IntegrateLayers(
  object = ic_slim,
  method = scVIIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.scvi",
  conda_env = "~/miniconda3/envs/scvi",
  group.by = "sample",
  verbose = TRUE
)

ic_slim <- RunUMAP(ic_slim, reduction = "harmony", reduction.name = "umap", dims = 1:30)
ic_slim <- RunUMAP(ic_slim, reduction = "rpca", reduction.name = "umap.rpca", dims = 1:30)
ic_slim <- RunUMAP(ic_slim, reduction = "integrated.scvi", reduction.name = "umap.scvi", dims = 1:30)

qsave(ic_slim, file.path("objects", "ic_reclustered.qs"))

# plot umaps ----
# rpca looked best
umap_ic_map <-
  DimPlot(ic, reduction = "umap.rpca", pt.size = .5, raster = FALSE, alpha = 0.1, group.by = "ic_level2", cols = my_cols_25, label = TRUE) +
  theme_rect() +
  xlab("UMAP1") +
  ylab("UMAP2")

ggsave(plot = umap_ic_map, file.path("results", "map", "ic_azimuth_rpca.png"), width = 10, height = 8)

Idents(ic_slim) <- ic_slim$ic_level2
umap_ic_map_split <-
  DimPlot(ic_slim, reduction = "umap.rpca", pt.size = .5, raster = FALSE, alpha = 0.1, split.by = "ic_level2", cols = my_cols_25, label = TRUE) +
  theme_rect()

ggsave(plot = umap_ic_map_split, file.path("results", "map", "ic_azimuth_rpca_split.png"), width = 40, height = 8)

# Optimize clustering resolution ----
ic_slim <- FindNeighbors(ic_slim, reduction = "rpca", dims = 1:30, assay = "RNA")

for (res in seq(from = 0.2, to = 2.3, by = 0.1)) {
  ic_slim <- FindClusters(ic_slim, resolution = res)
}

str(ic_slim@meta.data)

# plot clustering ---
resolutions <- paste0("RNA_snn_res.", seq(from = 0.2, to = 2.3, by = 0.1))

umap_list <- list()
for (res in resolutions) {
  umap_list[[res]] <-
    DimPlot(ic_slim, reduction = "umap.rpca", pt.size = .1, raster = FALSE, alpha = 0.1, group.by = res, cols = my_cols_50, label = TRUE) +
    theme_rect() +
    NoLegend()
}

umap_list <- patchwork::wrap_plots(umap_list, ncol = 2)
ggsave(plot = umap_list, file.path("results", "umap", "ic_rpca_umap_resolutions.png"), width = 16, height = 40)

umap_ic <-
  DimPlot(ic_slim, reduction = "umap.rpca", pt.size = .1, raster = FALSE, alpha = 0.1, group.by = "RNA_snn_res.2.3", cols = my_cols_50, label = TRUE) +
  theme_rect()

ggsave(plot = umap_ic, file.path("results", "umap", "ic_rpca_umap.png"), width = 10, height = 8)

library(clustree)
clustree <- clustree(ic_slim, prefix = "RNA_snn_res.")
ggsave(plot = clustree, file.path("results", "umap", "ic_clustree.png"), width = 15, height = 15)

# Visualize key marker genes with feature plots ----
endo_epi_macs <-
  FeaturePlot(
    ic,
    features = c("MS4A7", "CX3CR1", "TREM2", "LYVE1", "FOLR2", "TIMD4"),
    reduction = "umap.rpca",
    pt.size = 0.1,
    raster = FALSE,
    # alpha = 0.2,
    coord.fixed = TRUE,
    cols = c("#F0F0F0", "#CB181D"),
    order = TRUE,
    ncol = 3
  ) &
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1.0)
    ) &
    xlab("UMAP1") &
    ylab("UMAP2")
ggsave(file.path("results", "featureplot", "ic_endo_epi_macs.png"), endo_epi_macs, width = 8, height = 6)
ggsave(file.path("results", "featureplot", "ic_endo_epi_macs.pdf"), endo_epi_macs, width = 8, height = 6)

b_igh <-
  FeaturePlot(
    ic,
    features = c("IGHM", "IGHD", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2"),
    reduction = "umap.rpca",
    pt.size = 0.1,
    raster = FALSE,
    coord.fixed = TRUE,
    cols = c("#F0F0F0", "#CB181D"),
    order = TRUE,
    ncol = 2
  ) &
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1.0)
    ) &
    xlab("UMAP1") &
    ylab("UMAP2")
ggsave(file.path("results", "featureplot", "ic_bc_igh.png"), b_igh, width = 8, height = 16)

# dotplots ---
b_plasma <- subset(ic, subset = ic_cluster %in% c("Plasma", "B"))

dp_b_plasma <-
  DotPlot(b_plasma,
    features = c("IGHM", "IGHD", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2"),
    scale = FALSE
  ) +
  viridis::scale_color_viridis(option = "viridis") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic")) +
  xlab("") +
  ylab("")

ggsave(plot = dp_b_plasma, file.path("results", "dotplot", "dp_ic_bc_igh.pdf"), width = 4.5, height = 1.5)

dotPlot(
  path = file.path("lookup", "markers.csv"),
  object = b_plasma,
  par = "Bc_Ig",
  dot_min = 0.01,
  height = 2,
  width = 6
)

dotPlot(
  path = file.path("lookup", "markers.csv"),
  object = ic,
  par = "Bc_Ig",
  dot_min = 0.01,
  height = 8,
  width = 6
)

dotPlot(
  path = file.path("lookup", "markers.csv"),
  object = ic,
  par = "dotplot_ic_short",
  dot_min = 0.01,
  height = 7.5,
  width = 12
)

dotPlot(
  path = file.path("lookup", "markers.csv"),
  object = ic,
  par = "Macro18",
  dot_min = 0.01,
  height = 7.5,
  width = 20
)

# find markers ---
ic_slim <- JoinLayers(ic_slim)

# find markers helper function
findMarkers <- function(ident1, ident2 = NULL, object, only_pos, min_pct, logfc_threshold, assay = assay) {
  result <- Seurat::FindMarkers(object, ident.1 = ident1, ident.2 = ident2, min.pct = min_pct, logfc.threshold = logfc_threshold, only.pos = only_pos, assay = assay) |>
    tibble::rownames_to_column("gene") |>
    dplyr::filter(p_val_adj < 0.05) |>
    dplyr::relocate(gene, avg_log2FC, p_val, p_val_adj) |>
    dplyr::arrange(desc(avg_log2FC))
  return(result)
}

topmarkers_ic <-
  lapply(
    ic@misc$ic_cluster_order,
    function(x) {
      message("Processing cluster ", x)
      try(findMarkers(ident1 = x, object = ic, only_pos = TRUE, min_pct = 0.1, logfc_threshold = 0.25, assay = "RNA"))
    }
  )

names(topmarkers_ic) <- ic@misc$ic_cluster_order

writexl::write_xlsx(topmarkers_ic, file.path("results", "de", "topmarkers_ic.xlsx"))

# Additional reference mapping ----
# map to tabula sapiens blood ----
ts_blood <- qread("/home/mischko/Documents/beruf/forschung/scRNA_reference/tabula_sapiens/ts_blood.qs", nthread = 4)
ts_bm <- qread("/home/mischko/Documents/beruf/forschung/scRNA_reference/tabula_sapiens/ts_bm.qs", nthread = 4)

# function to map project query on ref and make predictions based on Seurat integration
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

predictions_prev <- mapSeurat(ref = ts_blood, query = ic_slim)
predictions_prev <- mapSeurat(ref = ts_bm, query = ic_slim)

# function to store predictions in seurat object
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

ic_slim <- storePred(predictions_prev, label_col = "ts_blood_label", score_col = "ts_blood_score", seu_obj = ic_slim)
ic_slim <- storePred(predictions_prev, label_col = "ts_bm_label", score_col = "ts_bm_score", seu_obj = ic_slim)
str(ic_slim@meta.data)

pred_plot_ts_blood <-
  DimPlot(ic_slim, reduction = "umap.rpca", group.by = "ts_blood_label", raster = FALSE, pt.size = .1, alpha = .1, cols = my_cols_25, label = TRUE) +
  theme_rect()
ggsave(plot = pred_plot_ts_blood, file.path("results", "map", "map_ic_ts_blood.png"), width = 20, height = 7)

pred_plot_ts_bm <-
  DimPlot(ic_slim, reduction = "umap.rpca", group.by = "ts_bm_label", raster = FALSE, pt.size = .1, alpha = .1, cols = my_cols_25, label = TRUE) +
  theme_rect()
ggsave(plot = pred_plot_ts_bm, file.path("results", "map", "map_ic_ts_bm.png"), width = 20, height = 7)

# run Azimuth bone marrow ---
options(timeout = 600)
ic_slim_bm <- RunAzimuth(query = ic_slim, reference = "bonemarrowref", assay = "RNA")

# prepare predictions
predictions_bm <-
  data.frame(
    bm_level1 = ic_slim_bm$predicted.celltype.l1,
    bm_level1_score = ic_slim_bm$predicted.celltype.l1.score,
    bm_level2 = ic_slim_bm$predicted.celltype.l2,
    bm_level2_score = ic_slim_bm$predicted.celltype.l2.score
  ) |>
  mutate(bm_level1 = ifelse(bm_level1_score < 0.4, "unknown", bm_level1)) |>
  mutate(bm_level2 = ifelse(bm_level2_score < 0.4, "unknown", bm_level2))

# add to seurat object
ic_slim <- AddMetaData(ic_slim, predictions_bm)

# save object
qsave(ic_slim, file.path("objects", "ic_slim.qs"))
qsave(ic_slim_bm, file.path("objects", "ic_slim_bm.qs"))

# plot predictions
pred_plot_azimuth_bm_level2 <-
  DimPlot(ic_slim, reduction = "umap.rpca", group.by = "bm_level2", raster = FALSE, pt.size = .1, alpha = .1, cols = my_cols_25, label = TRUE) +
  theme_rect()
ggsave(plot = pred_plot_azimuth_bm_level2, file.path("results", "map", "map_ic_azimuth_bm_level2.png"), width = 20, height = 7)

# run Azimuth tonsil ---
ic_slim_tonsil <- RunAzimuth(query = ic_slim, reference = "tonsilref", assay = "RNA")

# prepare predictions
predictions_tonsil <-
  data.frame(
    tonsil_level1 = ic_slim_tonsil$predicted.celltype.l1,
    tonsil_level1_score = ic_slim_tonsil$predicted.celltype.l1.score,
    tonsil_level2 = ic_slim_tonsil$predicted.celltype.l2,
    tonsil_level2_score = ic_slim_tonsil$predicted.celltype.l2.score
  ) |>
  mutate(tonsil_level1 = ifelse(tonsil_level1_score < 0.4, "unknown", tonsil_level1)) |>
  mutate(tonsil_level2 = ifelse(tonsil_level2_score < 0.4, "unknown", tonsil_level2))

# add to seurat object
ic_slim <- AddMetaData(ic_slim, predictions_tonsil)

# save object
qsave(ic_slim, file.path("objects", "ic_slim.qs"))
qsave(ic_slim_tonsil, file.path("objects", "ic_slim_tonsil.qs"))

# plot predictions
pred_plot_azimuth_tonsil_level1 <-
  DimPlot(ic_slim, reduction = "umap.rpca", group.by = "tonsil_level1", raster = FALSE, pt.size = .1, alpha = .1, cols = my_cols_25, label = TRUE) +
  theme_rect()
ggsave(plot = pred_plot_azimuth_tonsil_level1, file.path("results", "map", "map_ic_azimuth_tonsil_level1.png"), width = 20, height = 7)

pred_plot_azimuth_tonsil_level2 <-
  DimPlot(ic_slim, reduction = "umap.rpca", group.by = "tonsil_level2", raster = FALSE, pt.size = .1, alpha = .1, cols = my_cols_50, label = TRUE) +
  theme_rect()
ggsave(plot = pred_plot_azimuth_tonsil_level2, file.path("results", "map", "map_ic_azimuth_tonsil_level2.png"), width = 20, height = 7)

# run Azimuth lung ---
ic_slim_lung <- RunAzimuth(query = ic_slim, reference = "lungref", assay = "RNA")

str(ic_slim_lung@meta.data)

# prepare predictions
predictions_lung <-
  data.frame(
    lung_level1 = ic_slim_lung$predicted.ann_level_1,
    lung_level1_score = ic_slim_lung$predicted.ann_level_1.score,
    lung_level2 = ic_slim_lung$predicted.ann_level_2,
    lung_level2_score = ic_slim_lung$predicted.ann_level_2.score,
    lung_level3 = ic_slim_lung$predicted.ann_level_3,
    lung_level3_score = ic_slim_lung$predicted.ann_level_3.score
  ) |>
  mutate(lung_level1 = ifelse(lung_level1_score < 0.4, "unknown", lung_level1)) |>
  mutate(lung_level2 = ifelse(lung_level2_score < 0.4, "unknown", lung_level2)) |>
  mutate(lung_level3 = ifelse(lung_level3_score < 0.4, "unknown", lung_level3))

# add to seurat object
ic_slim <- AddMetaData(ic_slim, predictions_lung)

# save object
qsave(ic_slim, file.path("objects", "ic_slim.qs"))
qsave(ic_slim_lung, file.path("objects", "ic_slim_lung.qs"))

# plot predictions
pred_plot_azimuth_lung_level3 <-
  DimPlot(ic_slim, reduction = "umap.rpca", group.by = "lung_level3", raster = FALSE, pt.size = .1, alpha = .1, cols = my_cols_25, label = TRUE) +
  theme_rect()
ggsave(plot = pred_plot_azimuth_lung_level3, file.path("results", "map", "map_ic_azimuth_lung_level3.png"), width = 20, height = 7)

# Final cluster annotation ----
# annotate immune cell clusters ----
# Manual annotation based on marker genes and reference mapping
Idents(ic_slim) <- ic_slim$RNA_snn_res.2.3

annotations <- readxl::read_xlsx(file.path("lookup", "ic_cluster_annotation.xlsx")) |>
  mutate(cluster = as.character(cluster))

names(annotations$final) <- levels(ic_slim)
ic_slim <- RenameIdents(ic_slim, annotations$final)

# save custom cluster and condition order and colors in seurat object
ic_slim@misc$ic_cluster_order <- annotations$cluster_order[!is.na(annotations$cluster_order)]
ic_slim@misc$ic_cluster_col <- setNames(my_cols_50[1:length(ic_slim@misc$ic_cluster_order)], ic_slim@misc$ic_cluster_order)

# sanity check
setdiff(ic_slim@misc$ic_cluster_order, levels(ic_slim))
setdiff(levels(ic_slim), ic_slim@misc$ic_cluster_order)
any(duplicated(ic_slim@misc$ic_cluster_order))

# save ordered annotations in meta.data
ic <- subset(ic_slim, idents = "remove", invert = TRUE)

ic$ic_cluster <- factor(Idents(ic), levels = ic@misc$ic_cluster_order)
Idents(ic) <- ic$ic_cluster

umap_ic <-
  DimPlot(ic, reduction = "umap.rpca", pt.size = .1, raster = FALSE, alpha = 0.3, group.by = "ic_cluster", cols = ic@misc$ic_cluster_col, label = TRUE) +
  theme_rect() +
  NoLegend() +
  xlab("UMAP1") +
  ylab("UMAP2")

ggsave(plot = umap_ic, file.path("results", "umap", "ic_rpca_annotated.png"), width = 8, height = 7)

# Save final annotated object ----
qs::qsave(ic, file.path("objects", "ic.qs"))

# compare macro18 markers with samc ---
macro18_markers <- findMarkers(ident1 = "Macro18", object = ic, only_pos = TRUE, min_pct = 0.1, logfc_threshold = 0.25, assay = "RNA")

# load data from stroke project
stroke <- qs::qread(file.path("objects", "stroke.qs"))

# need to normalize RNA assay
DefaultAssay(stroke) <- "RNA"
stroke <- NormalizeData(stroke)

# use same settings for samc markers
samc_markers <- findMarkers(ident1 = "SAMC", object = stroke, only_pos = TRUE, min_pct = 0.1, logfc_threshold = 0.25, assay = "RNA")

lookup_samc_markers <- homologene::mouse2human(samc_markers$gene)

samc_markers_ortho <-
  samc_markers |>
  left_join(lookup_samc_markers, by = c("gene" = "mouseGene"))

samc_top <-
  samc_markers_ortho |>
  slice(1:100) |>
  dplyr::filter(humanGene %in% rownames(ic)) |>
  pull(humanGene) |>
  na.omit()

ic <- AddModuleScore(
  object = ic,
  features = list(samc_top),
  name = "SAMC"
)

qs::qsave(ic, file.path("objects", "ic.qs"))

# better name
ic$SAMC <- ic$SAMC1
ic$SAMC1 <- NULL

fplot_samc <-
  FeaturePlot(
    ic,
    features = "SAMC",
    reduction = "umap.rpca",
    pt.size = 0.1,
    raster = FALSE,
    label = TRUE,
    coord.fixed = TRUE,
    order = TRUE
  ) +
    viridis::scale_color_viridis() +
    theme_rect() +
    xlab("UMAP1") &
    ylab("UMAP2")

ggsave(file.path("results", "module", "ic_samc.png"), fplot_samc, width = 8, height = 6)

# map SAMC from stroke to immune cells
# use functions from  integrate.R
stroke_human <- convertRownames(stroke)
predictions_stroke <- mapSeurat(ref = stroke_human, query = ic)
ic <- storePred(predictions_stroke, label_col = "stroke_label", score_col = "stroke_score", seu_obj = ic)

# plot predicted SAMC
Idents(ic) <- ic$stroke_label

# Function to plot predicted SAMC
plot_predicted_samc <- function(seu_obj, output_path) {
  # Create a data frame with alpha values, color values, and UMAP coordinates
  umap_coords <- Embeddings(seu_obj, "umap.rpca")
  alpha_values <- ifelse(seu_obj$stroke_label == "SAMC", 0.5, 0.001)
  color_values <- ifelse(seu_obj$stroke_label == "SAMC", "SAMC", "Other")
  color_map <- setNames(c("red", "grey"), c("SAMC", "Other"))
  plot_df <- data.frame(
    cell = rownames(umap_coords),
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    alpha = alpha_values,
    color = color_values
  )

  # Create custom ggplot
  samc_predicted <- ggplot(
    plot_df,
    aes(x = UMAP1, y = UMAP2, alpha = alpha, color = color)
  ) +
    geom_point(size = 0.1) +
    scale_color_manual(values = color_map) + # Simplified color scale
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1.0),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "none"
    ) +
    ggtitle("Predicted SAMC")

  ggsave(output_path, samc_predicted, width = 4, height = 4)
}

# plot predicted SAMC
Idents(ic) <- ic$stroke_label
plot_predicted_samc(ic, file.path("results", "map", "samc_predicted.png"))

Idents(ic) <- ic$ic_cluster
qsave(ic, file.path("objects", "ic.qs"))

# featureplot of SAMC markers
ic_samc_markers_fplot <-
  FeaturePlot(
    ic,
    features = c("SPP1", "APOE", "LPL", "FABP5", "GPNMB"),
    reduction = "umap.rpca",
    cols = c("#F0F0F0", "#CB181D"),
    pt.size = 0.1,
    raster = FALSE,
    coord.fixed = TRUE,
    order = TRUE
  ) &
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1.0)
    ) &
    xlab("UMAP1") &
    ylab("UMAP2")

ggsave(file.path("results", "featureplot", "ic_samc_markers.png"), ic_samc_markers_fplot, width = 8, height = 12)


# dotplot of SAMC markers
Idents(ic) <- ic$ic_cluster
dotPlot(
  path = file.path("lookup", "markers.csv"),
  object = ic,
  par = "SAMC",
  dot_min = 0.01,
  height = 6,
  width = 4.5
)

