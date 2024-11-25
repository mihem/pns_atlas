#===============================================================================
# Cerebro Export Script
#===============================================================================
# Purpose: Prepare and export Seurat object for CerebroApp visualization:
# - Convert Seurat object to Cerebro-compatible format
# - Clean up unnecessary data and optimize memory usage
# - Export selected metadata and marker genes
# - Generate interactive visualization file
#===============================================================================

# Install custom cerebroApp version ----
remotes::install_github("mihem/cerebroAppLite")

# load libraries ----
library(Seurat)
library(BPCells)
library(SeuratObject)
library(tidyverse)
library(conflicted)
library(qs)
library(cerebroAppLite)

# general settings ----
options(warn = 0) 
conflicts_prefer(base::setdiff)

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)

# prepare data for cerebro ----
# Downsample and clean up Seurat object
sc_merge_small <- subset(sc_merge, downsample = 1000)

sc_merge_cerebro <- sc_merge

# Remove unnecessary data to optimize memory usage
sc_merge_cerebro <- DietSeurat(
    sc_merge_cerebro,
    layers = c("data"),
    assays = "RNA",
    dimreducs = c("umap.scvi.full")
)

# Additional cleanup
sc_merge_cerebro[["RNA"]]$scale.data <- NULL
sc_merge_cerebro[["RNA"]]$counts <- NULL

# # convert to in memory matrix
# sc_merge_cerebro[["RNA"]]$counts <- as(object = sc_merge_cerebro[["RNA"]]$counts, Class = "dgCMatrix")

# # convert to v3 assay
# sc_merge_cerebro[["RNA"]] <- as(object = sc_merge_cerebro[["RNA"]], Class = "Assay")


# Prepare metadata for export ----
# Select and rename relevant metadata columns
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

# Calculate marker genes ----
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

# Export for CerebroApp ----
# Save processed object and export to Cerebro format
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

# Launch CerebroApp ----
cerebroAppLite::launchCerebro()