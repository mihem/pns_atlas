# prepare data for cerebroApp

# install my version of cerebroApp
remotes::install_github("mihem/cerebroAppLite")

# load libraries
library(Seurat)
library(BPCells)
library(SeuratObject)
library(tidyverse)
library(conflicted)
library(qs)
library(cerebroAppLite)

# general settings  ----
options(warn = 0)
conflicts_prefer(base::setdiff)

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)

# downsample ----
sc_merge_small <- subset(sc_merge, downsample = 1000)

##export from seurat to cerebro format
sc_merge_cerebro <- sc_merge
sc_merge_cerebro <- sc_merge_small

# remove unnecessary assays and dims
sc_merge_cerebro <- DietSeurat(
    sc_merge_cerebro,
    layers = c("data"),
    assays = "RNA",
    dimreducs = c("umap.scvi.full")
)

str(sc_merge_cerebro$RNA$counts@matrix, max.level = 3)

# somehow the scale.data is still present, only data required for gene expression
sc_merge_cerebro[["RNA"]]$scale.data <- NULL
sc_merge_cerebro[["RNA"]]$counts <- NULL

# write_matrix_dir(sc_merge_cerebro$RNA$data, dir = paste0("matrix_cerebro/"))

# # convert to in memory matrix
# sc_merge_cerebro[["RNA"]]$counts <- as(object = sc_merge_cerebro[["RNA"]]$counts, Class = "dgCMatrix")

# # convert to v3 assay
# sc_merge_cerebro[["RNA"]] <- as(object = sc_merge_cerebro[["RNA"]], Class = "Assay")

# only keep relevant metadata
sc_merge_cerebro@meta.data <-
    sc_merge_cerebro@meta.data |>
    tibble::rownames_to_column("barcode") |>
    rename(
        diagnosis = level2,
        group = level0,
        suter_p60 = suter_p60_label_full,
        milbrandt_sciatic = milbrandt_sciatic_label_full
    ) |>
    select(barcode, cluster, group, diagnosis, nCount_RNA, nFeature_RNA, percent_mt, sample, sex, age, center, milbrandt_sciatic, suter_p60) |>
    tibble::column_to_rownames(var = "barcode")

sc_merge_cerebro <- getMarkerGenes(
  sc_merge_cerebro,
  organism = "hg",
  groups = c("cluster"),
  only_pos = TRUE,
  assay = "RNA",
  min_pct = 0.1,
  thres_logFC = 0.25,
  thres_p_val = 0.05
)

qs::qsave(sc_merge_cerebro, file = file.path("objects", "sc_merge_cerebro.qs"))

cerebroAppLite::exportFromSeurat(
  object = sc_merge_cerebro,
  file = file.path("objects", "sc_merge_cerebro.crb"),
  experiment_name = "Sural",
  groups = c("cluster", "sample", "group", "diagnosis"),
  organism = "hg",
  nUMI = "nCount_RNA",
  nGene = "nFeature_RNA",
  use_delayed_array = FALSE
)

cerebroAppLite::launchCerebro()