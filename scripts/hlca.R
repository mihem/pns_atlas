# libraries  ----
library(Seurat)
library(BPCells)
library(SeuratObject)
library(tidyverse)
library(patchwork)
library(conflicted)
library(qs)
library(pals)
library(scMisc)
library(Polychrome)
library(Azimuth)
library(SeuratData)

# azimuth map to human lung cell atlas ----
available_data <- SeuratData::AvailableData()

# higher timeout
options(timeout = 600)

# DefaultAssay(sc_merge) <- "sketch"
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)
sc_small <- qs::qread(file.path("objects", "sc_small.qs"), nthread = 4)
sc_small$RNA$counts <- as(object = sc_small$RNA$counts, Class = "dgCMatrix")
sc_small$RNA$data <- as(object = sc_small$RNA$data, Class = "dgCMatrix")

sc_small <- subset(sc_merge, downsample = 1000)

str(sc_small$RNA@layers, max.level = 1)
str(sc_small$sketch@layers, max.level = 1)

# DefaultAssay(sc_merge) <- "RNA"
# str(sc_merge_matrix@assays$RNA@layers, max.level = 2)
# str(sc_merge_matrix@assays$sketch@layers, max.level = 2)

# sc_merge_matrix$RNA$counts <- as(object = sc_merge_matrix$RNA$counts, Class = "dgCMatrix")
# sc_merge_matrix$RNA$data <- as(object = sc_merge_matrix$RNA$data, Class = "dgCMatrix")
# sc_merge_matrix <- ScaleData(sc_merge_matrix)

# sc_merge_ec <- RunAzimuth(query = sc_merge_matrix, reference = "lungref", assay = "sketch")

sc_small_hlca <- RunAzimuth(query = sc_small, reference = "lungref", assay = "sketch")
qs::qsave(sc_small_hlca, file.path("objects", "sc_small_hlca.qs"))

DefaultAssay(sc_merge) <- "sketch"

DefaultAssay(sc_merge) <- "RNA"
sc_merge <- split(x = sc_merge, f = sc_merge$sample)


for (i in Layers(sc_merge)){
    sc_merge[["RNA"]][[i]] <- as(sc_merge[["RNA"]][[i]], "dgCMatrix")
    sc_merge[["RNA"]][[i]] <- write_matrix_dir(sc_merge[["RNA"]][[i]], dir = paste0("matrix_final2/", gsub("counts.", "", i)))
}

qs::qsave(sc_merge, file.path("objects", "sc_merge_2.qs"))

sc_merge <- JoinLayers(sc_merge)
sc_hlca <- RunAzimuth(query = sc_merge, reference = "lungref")

str(sc_merge@assays, max.level = 4)

str(sc_merge$RNA$counts@matrix@matrix@matrix_list[[1]]@matrix@matrix_list[[1]]@matrix@matrix, max.level = 2)
str(sc_merge$RNA$counts@matrix@matrix, max.level = 2)
str(sc_merge$RNA$counts@matrix@matrix@matrix_list, max.level = 3)

# fp1 <- FeaturePlot(sc_small_ec, features = "predictionscoreannlevel3_Macrophages", reduction = "umap.scvi.full")
# ggsave(plot = fp1, "test3.png")

# fp2 <- FeaturePlot(sc_small_ec, features = "predictionscoreannlevel3_EC venous", reduction = "umap.scvi.full")
# ggsave(plot = fp2, "test4.png")

# prepare prediction 
predictions_hlca <-
  data.frame(
    hlca_lev3 = sc_small_hlca$predicted.ann_level_3,
    hlca_lev3_score = sc_small_hlca$predicted.ann_level_3.score,
    hlca_lev4 = sc_small_hlca$predicted.ann_level_4,
    hlca_lev4_score = sc_small_hlca$predicted.ann_level_4.score
  ) |>
  mutate(hlca_lev3 = ifelse(hlca_lev3_score < 0.4, "unknown", hlca_lev3)) |>
  mutate(hlca_lev4 = ifelse(hlca_lev4_score < 0.4, "unknown", hlca_lev4))

# add to seurat object
sc_small <- AddMetaData(sc_small, predictions_hlca)

# plot prediction from HLCA
pred_plot_hlca_lev3_small <-
 DimPlot(sc_small, reduction = "umap.scvi.full", group.by = "hlca_lev3", raster = FALSE, pt.size = .1, alpha = .1, cols = rev(my_cols_25), label = TRUE) +
  theme_rect()
ggsave(plot = pred_plot_hlca_lev3_small, file.path("results", "map", "map_hlca_lev3_small.png"), width = 15, height = 10)

qs::qsave(sc_small, file.path("objects", "sc_small.qs"))

renv::hydrate("cellxgene.census")
library(cellxgene.census)

census <-  open_soma()

organism <-  "Homo sapiens"
gene_filter <-  "feature_id %in% c('ENSG00000107317', 'ENSG00000106034')"
cell_filter <-   "cell_type == 'sympathetic neuron'"
cell_columns <-  c("assay", "cell_type", "tissue", "tissue_general", "suspension_type", "disease")

seurat_obj <-  get_seurat(
   census = census,
   organism = organism,
   var_value_filter = gene_filter,
   obs_value_filter = cell_filter,
   obs_column_names = cell_columns
)

Error in match.arg(arg = layer, choices = Layers(object = object, search = FALSE)) : 
  'arg' should be one of “counts”, “data”, “scale.data”

9: stop(sprintf(ngettext(length(chs <- unique(choices[nzchar(choices)])), 
       "'arg' should be %s", "'arg' should be one of %s"), paste(dQuote(chs), 
       collapse = ", ")), domain = NA)
8: match.arg(arg = layer, choices = Layers(object = object, search = FALSE))
7: `LayerData<-.Assay`(object = `*tmp*`, layer = i, ..., value = value)
6: `LayerData<-`(object = `*tmp*`, layer = i, ..., value = value)
5: `[[<-`(`*tmp*`, names(var), value = structure(list(feature_name = c("CPED1", 
   "PTGDS"), feature_length = c(7683L, 2712L)), row.names = c("ENSG00000106034", 
   "ENSG00000107317"), class = "data.frame"))
4: `[[<-`(`*tmp*`, names(var), value = structure(list(feature_name = c("CPED1", 
   "PTGDS"), feature_length = c(7683L, 2712L)), row.names = c("ENSG00000106034", 
   "ENSG00000107317"), class = "data.frame"))
3: self$to_seurat_assay(X_layers = X_layers, obs_index = obs_index, 
       var_index = var_index, var_column_names = var_column_names)
2: expt_query$to_seurat(X_layers = X_layers, obs_column_names = obs_column_names, 
       var_column_names = var_column_names, var_index = var_index)
1: get_seurat(census = census, organism = organism, var_value_filter = gene_filter, 
       obs_value_filter = cell_filter, obs_column_names = cell_columns)

# Open obs SOMADataFrame
cell_metadata <-  census$get("census_data")$get("homo_sapiens")$get("obs")

# Read as Arrow Table
cell_metadata <-  cell_metadata$read(
   value_filter = "sex == 'female' & cell_type %in% c('microglial cell', 'neuron')",
   column_names = c("assay", "cell_type", "sex", "tissue", "tissue_general", "suspension_type", "disease")
)

# Concatenates results to an Arrow Table
cell_metadata <-  cell_metadata$concat()

# Convert to R tibble (dataframe)
cell_metadata <-  as.data.frame(cell_metadata)

reference <- LoadData("lungref", type = "azimuth")$map

anchors <- FindTransferAnchors(
  reference = reference,
  query = sc_merge,
  k.filter = NA,
  reference.neighbors = "refdr.annoy.neighbors",
  reference.assay = "refAssay",
  query.assay = "sketch",
  reference.reduction = "refDR",
  normalization.method = "SCT",
  features = rownames(Loadings(reference[["refDR"]]))
)


data("pbmc3k")
reference <- UpdateSeuratObject(pbmc3k)
reference_list <- list(ref = reference, query = sc_merge)
features <- SelectIntegrationFeatures(object.list = reference_list)

anchors <- FindTransferAnchors(
    reference = reference_list$ref,
    query = reference_list$query,
    features = features,
    query.assay = "RNA"
  )

predictions <- TransferData(anchorset = anchors, refdata = pbmc3k$seurat_annotations)

sc_merge <- storePred(predictions, label_col = "pbmc3k_label_2", score_col = "pbmc3k_score_2", seu_obj = sc_merge)
dplyr::count(sc_merge@meta.data, pbmc_3k_label)
dplyr::count(sc_merge@meta.data, pbmc3k_label_2)

pred_plot_pbmc3k_sketch <-
  DimPlot(sc_merge, reduction = "umap.scvi.full", group.by = "pbmc_3k_label", raster = FALSE, pt.size = .1, alpha = .1, cols = my_cols_25, label = TRUE) +
  theme_rect()
ggsave(plot = pred_plot_pbmc3k_sketch, file.path("results", "map", "map_pbmc3k_sketch.png"), width = 8, height = 8)

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

sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthreads = 4)
sc_merge[["sketch"]] <- JoinLayers(sc_merge[["sketch"]])

DefaultAssay(sc_merge) <- "sketch"

predictions_prev <- mapSeurat(ref = sc_merge_prev_small, query = sc_merge)
predictions_milbrandt <- mapSeurat(ref = human_pns_sciatic_milbrandt, query = sc_merge)
