# subset immune cells and characterize them in more depth

# libraries  ----
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

# devtools::install_github("mihem/scMisc")
# detach(package:scMisc, unload = TRUE)

# general settings  ----
options(warn = 0)
options(Seurat.object.assay.version = "v5")
future::plan("multicore", workers = 6)
options(future.globals.maxSize = 16000 * 1024^2)
conflicts_prefer(base::setdiff)
my_cols_25 <- pals::cols25()
my_cols_50 <- unname(Polychrome::createPalette(50, pals::cols25()))

# load preprocessed data ----
ic <- qread(file.path("objects", "ic.qs"), nthread = 4)
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)

# subset immune cells ----
ic <- subset(sc_merge, cluster %in% c("Macro1", "Macro2", "Granulo", "B", "T_NK", "Mast"))

ic$cluster <- droplevels(ic$cluster)

DimPlot(ic, reduction = "umap.scvi.full", group.by = "cluster", raster = FALSE, pt.size = .1, alpha = .1, label = TRUE)

# run Azimuth ---
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
    )|>
  mutate(ic_level1 = ifelse(ic_level1_score < 0.4, "unknown", ic_level1)) |>
  mutate(ic_level2 = ifelse(ic_level2_score < 0.4, "unknown", ic_level2)) |>
  mutate(ic_level3 = ifelse(ic_level3_score < 0.4, "unknown", ic_level3))

# add to seurat object ---
ic <- AddMetaData(ic, predictions_pbmcref)

# plot predictions ----
pred_plot_ic_lev2 <-
 DimPlot(ic, reduction = "umap.scvi.full", group.by = "ic_level2", raster = FALSE, pt.size = .1, alpha = .1, cols = rev(my_cols_25), label = FALSE )+
  theme_rect()

ggsave(plot = pred_plot_ic_lev2, file.path("results", "map", "map_ic_lev2.png"), width = 15, height = 10)


# remove unnecessary layers and assays ----
ic_slim <- DietSeurat(
  ic,
  assay = c("RNA"),
  layers = c("counts"),
)

# convert to sparse matrices ----
ic_slim$RNA$counts <- as(object = ic_slim$RNA$counts, Class = "dgCMatrix")
ic_slim$RNA$data <- as(object = ic_slim$RNA$data, Class = "dgCMatrix")
ic_slim[["RNA"]]$scale.data <- NULL

# integrate ----
ic_slim[["RNA"]] <- split(
  x = ic_slim[["RNA"]],
  f = ic_slim$sample)

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
    DimPlot(ic_slim, reduction = "umap.rpca", pt.size = .5, raster = FALSE, alpha = 0.1, group.by = "ic_level2", cols = my_cols_25, label = TRUE) +
    theme_rect()

ggsave(plot = umap_ic_map, file.path("results", "map", "ic_azimuth_rpca.png"), width = 10, height = 8)

Idents(ic_slim) <- ic_slim$ic_level2
umap_ic_map_split <-
    DimPlot(ic_slim, reduction = "umap.rpca", pt.size = .5, raster = FALSE, alpha = 0.1, split.by = "ic_level2", cols = my_cols_25, label = TRUE) +
    theme_rect()

ggsave(plot = umap_ic_map_split, file.path("results", "map", "ic_azimuth_rpca_split.png"), width = 40, height = 8)

ic_slim <- FindNeighbors(ic_slim, reduction = "rpca", dims = 1:30, assay = "RNA")

for (res in seq(from = 2.2, to = 2.3, by = 0.1)) {
  ic_slim <- FindClusters(ic_slim, resolution = res)
}

str(ic_slim@meta.data)

# plot clustering ---
resolutions <- paste0("RNA_snn_res.", seq(from = 0.2, to = 2.0, by = 0.1))

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


# feature plots ----
scMisc::fPlot(ic_slim,
  path = file.path("lookup", "markers.csv"),
  par = "cellmarkers_aie",
  reduction = "umap.rpca",
)

scMisc::fPlot(ic,
  path = file.path("lookup", "markers.csv"),
  par = "NK",
  reduction = "umap.rpca",
)

scMisc::fPlot(ic_slim,
  path = file.path("lookup", "markers.csv"),
  par = "tcells_subcluster_aie",
  reduction = "umap.rpca",
)

scMisc::fPlot(ic,
  path = file.path("lookup", "markers.csv"),
  par = "MAIT",
  reduction = "umap.rpca",
)

scMisc::fPlot(ic,
  path = file.path("lookup", "markers.csv"),
  par = "iNKT",
  reduction = "umap.rpca",
)

# features plot endoMacro epiMacro
macro_names <- str_subset(ic@misc$ic_cluster_order, "Macro")
macro <- subset(ic, subset = ic_cluster %in% macro_names)

scMisc::fPlot(
  macro,
  path = file.path("lookup", "markers.csv"),
  par = "EpiMac_Ydens",
  reduction = "umap.rpca",
)

scMisc::fPlot(
  macro,
  path = file.path("lookup", "markers.csv"),
  par = "EndoMac_Ydens",
  reduction = "umap.rpca",
)

FeaturePlot(ic, features = markers$MAIT, reduction = "umap.rpca", order = TRUE)
ggsave(file.path("results", "featureplot", "fp_ic_MAIT.pdf"), width = 10, height = 8)

FeaturePlot(ic, features = markers$iNKT, reduction = "umap.rpca", order = TRUE)
ggsave(file.path("results", "featureplot", "fp_ic_iNKT.pdf"), width = 10, height = 8)

# dotplots ---
dotPlot(
  path = file.path("lookup", "markers.csv"),
  object = ic,
  par = "dotplot_ic_short",
  dot_min = 0.01,
  height = 7.5,
  width = 12
) 

# epiMacs endoMacs
macro_names <- str_subset(ic@misc$ic_cluster_order, "Macro")
macro <- subset(ic, subset = ic_cluster %in% macro_names)

dotPlot(
  path = file.path("lookup", "markers.csv"),
  object = macro,
  par = "EpiMac_Ydens",
  dot_min = 0.01,
  height = 6.5,
  width = 5.5,
) 

dotPlot(
  path = file.path("lookup", "markers.csv"),
  object = macro,
  par = "EndoMac_Ydens",
  dot_min = 0.01,
  height = 6.5,
  width = 8,
) 


lapply(
  paste0("DC", 1:6, "_villani"),
  FUN = function(x) {
    dotPlot(
      path = file.path("lookup", "markers.csv"),
      object = ic_slim,
      par = x,
      dot_min = 0.01,
      height = 8,
      width = 10
    )
  })

dotPlot(
  path = file.path("lookup", "markers.csv"),
  object = ic_slim,
  par = "dotplot_short",
  dot_min = 0.01,
  height = 8,
  width = 26
) 

dotPlot(
  path = file.path("lookup", "markers.csv"),
  object = ic_slim,
  par = "PNS_Macro",
  dot_min = 0.01,
  height = 8,
  width = 6
) 

dotPlot(
  path = file.path("lookup", "markers.csv"),
  object = ic_slim,
  par = "microglia",
  dot_min = 0.01,
  height = 8,
  width = 6
) 

# endo_epi_Macs_Ydens <- read_csv("/home/mischko/Documents/beruf/forschung/markers/markers_pns.csv") |>
#   dplyr::filter(cell_source %in% c("EndoMac_Ydens", "EpiMac_Ydens")) |>
#   dplyr::filter(!is.na(gene))

# print(endo_epi_Macs_Ydens, n = Inf)

lapply(
  c("EndoMac_Ydens", "EpiMac_Ydens"),
  FUN = function(x) {
    dotPlot(
      path = file.path("lookup", "markers.csv"),
      object = ic_slim,
      par = x,
      dot_min = 0.01,
      height = 8,
      width = 7
    )
  })

# join layers ---
ic_slim <- JoinLayers(ic_slim)

# find markers

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

# remove failed lists (because of too few cells)
# topmarkers_ic <- topmarkers_ic[sapply(topmarkers_ic, is.list)]

writexl::write_xlsx(topmarkers_ic, file.path("results", "de", "topmarkers_ic.xlsx"))

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

#function to store predictions in seurat object
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
  )   |>
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

# annotate immune cell clusters ----
Idents(ic_slim) <- ic_slim$RNA_snn_res.2.3

annotations <- readxl::read_xlsx(file.path("lookup", "ic_cluster_annotation.xlsx")) |>
    mutate(cluster = as.character(cluster)) 

names(annotations$final) <- levels(ic_slim)
ic_slim <- RenameIdents(ic_slim, annotations$final)

# save custom cluster and condition order and colors in seurat object
ic_slim@misc$ic_cluster_order <- annotations$cluster_order[!is.na(annotations$cluster_order)]
ic_slim@misc$ic_cluster_col <- setNames(my_cols_50[1:length(ic_slim@misc$ic_cluster_order)], ic_slim@misc$ic_cluster_order)

#sanity check
setdiff(ic_slim@misc$ic_cluster_order, levels(ic_slim))
setdiff(levels(ic_slim), ic_slim@misc$ic_cluster_order)
any(duplicated(ic_slim@misc$ic_cluster_order))

# save ordered annotations in meta.data
ic <- subset(ic_slim, idents = "remove", invert = TRUE)

ic$ic_cluster <- factor(Idents(ic), levels = ic@misc$ic_cluster_order)
Idents(ic) <- ic$ic_cluster

umap_ic <-
    DimPlot(ic, reduction = "umap.rpca", pt.size = .5, raster = FALSE, alpha = 0.1, group.by = "ic_cluster", cols = ic@misc$ic_cluster_col, label = TRUE) +
    theme_rect() +
    NoLegend()
ggsave(plot = umap_ic, file.path("results", "umap", "ic_rpca_annotated.png"), width = 12, height = 9)

# save object ---
qs::qsave(ic, file.path("objects", "ic.qs"))
