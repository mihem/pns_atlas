#===============================================================================
# Marker Gene Analysis and Cross-Dataset Comparison Script
#===============================================================================
# Purpose: Compare marker genes across multiple peripheral nerve datasets:
# - Our dataset (sc_merge)
# - Milbrandt sciatic nerve atlas
# - Suter P60 dataset
# - Carr mesenchymal dataset
# - Wolbert healthy dataset
#===============================================================================

# Libraries and Setup ----
library(Seurat)
library(BPCells)
library(SeuratObject)
library(tidyverse)
library(conflicted)
library(qs)
library(scMisc)
library(writexl)
library(readxl)
library(RColorBrewer)
library(homologene)

# general settings  ----
options(warn = 0)
future::plan("multicore", workers = 6)
conflicts_prefer(base::setdiff)

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)

# seurat map to reference datasets ----
pns_sn_sciatic_milbrandt <- qs::qread("/home/mischko/Documents/beruf/forschung/scRNA_reference/pns_atlas_milbrandt/pns_sn_sciatic_GSE182098.qs", nthreads = 4)
suter_p60 <- readRDS("/home/mischko/Documents/beruf/forschung/scRNA_reference/sciatic_nerve_atlas_suter/10xGeno_P60.rds")
carr_comMes <- qs::qread("/home/mischko/Documents/beruf/forschung/scRNA_reference/pns_carr/carr_combMes.qs")
load("/home/mischko/Documents/beruf/forschung/scRNA_reference/pns_carr/InjUninjMesenchymalCombined.RData")

wolbert <- qread("/home/mischko/Documents/beruf/forschung/scRNA_reference/pns_wolbert/wolbert_healthy.qs")

renv::install("BaderLab/scClustViz")

str(InjUninjMesenchymalCombined_sCVdL)

janitor::make_clean_names(names(topmarkers_suter))


#' Find Marker Genes for Cell Type Comparison
findMarkers <- function(ident1, ident2 = NULL, object, only_pos, min_pct, logfc_threshold, assay = assay) {
  result <- Seurat::FindMarkers(object, ident.1 = ident1, ident.2 = ident2, min.pct = min_pct, logfc.threshold = logfc_threshold, only.pos = only_pos, assay = assay) |>
    tibble::rownames_to_column("gene") |>
    dplyr::filter(p_val_adj < 0.05) |>
    dplyr::relocate(gene, avg_log2FC, p_val, p_val_adj) |>
    dplyr::arrange(desc(avg_log2FC))
  return(result)
}

# Marker Analysis Per Dataset ----
# Calculate markers for each cluster across datasets
topmarkers <-
  lapply(
    unique(sc_merge@misc$cluster_order),
    function(x) {
      message("Processing cluster ", x)
      try(findMarkers(ident1 = x, object = sc_merge, only_pos = TRUE, min_pct = 0.1, logfc_threshold = 0.25, assay = "RNA"))
    }
  )

# or read in previously calculated markers for our dataset
topmarkers <-
    lapply(
        sc_merge@misc$cluster_order,
        function(x) {
            read_excel(file.path("results", "de", "topmarkers_final.xlsx"), sheet = x)
        }
    ) |>
    setNames(sc_merge@misc$cluster_order)

# Reference Dataset Markers ----
# Calculate markers for Milbrandt, Suter, and Wolbert datasets
topmarkers_milbrandt <-
    lapply(
        unique(pns_sn_sciatic_milbrandt$cluster),
        function(x) {
            message("Processing cluster ", x)
            try(findMarkers(ident1 = x, object = pns_sn_sciatic_milbrandt, only_pos = TRUE, min_pct = 0.1, logfc_threshold = 0.25, assay = "RNA"))
        }
    ) |>
    setNames(unique(pns_sn_sciatic_milbrandt$cluster))

write_xlsx(topmarkers_milbrandt, file.path("results", "de", "topmarkers_milbrandt.xlsx"))

topmarkers_suter <-
    lapply(
        unique(suter_p60$cluster),
        function(x) {
            message("Processing cluster ", x)
            try(findMarkers(ident1 = x, object = suter_p60, only_pos = TRUE, min_pct = 0.1, logfc_threshold = 0.25, assay = "RNA"))
        }
    ) |>
    setNames(unique(suter_p60$cluster))

names(topmarkers_suter)[8] <- "Per_VSMC"
names(topmarkers_suter)[9] <- "Per_EC"

write_xlsx(topmarkers_suter, file.path("results", "de", "topmarkers_suter.xlsx"))

topmarkers_wolbert <-
    lapply(
        unique(wolbert$cluster),
        function(x) {
            message("Processing cluster ", x)
            try(findMarkers(ident1 = x, object = wolbert, only_pos = TRUE, min_pct = 0.1, logfc_threshold = 0.25, assay = "RNA"))
        }
    ) |>
    setNames(unique(wolbert$cluster))

write_xlsx(topmarkers_wolbert, file.path("results", "de", "topmarkers_wolbert.xlsx"))

# Reference Integration ----
# Make cluster names consistent and merge reference data
names(topmarkers_milbrandt)[names(topmarkers_milbrandt) == "EnFibro"] <- "endoC"
names(topmarkers_milbrandt)[names(topmarkers_milbrandt) == "EpC"] <- "epiC"
names(topmarkers_milbrandt)[names(topmarkers_milbrandt) == "PnC"] <- "periC"

Idents(sc_merge) <- sc_merge$cluster
topmarkers[["periC"]] <- findMarkers(
    ident1 = "periC",
    object = sc_merge,
    only_pos = TRUE,
    min_pct = 0.1,
    logfc_threshold = 0.25,
    assay = "RNA"
)

# Create comparative visualizations ----
sc_merge_subset <- subset(sc_merge, subset = cluster %in% c("mySC", "nmSC", "periC1", "periC2", "periC3"))

# rename periC1, periC2, periC3 to periC
sc_merge_subset$cluster <- gsub(pattern = "periC\\d", replacement = "periC", x = sc_merge_subset$cluster)
sc_merge_subset$cluster <- factor(sc_merge_subset$cluster, levels = c("mySC", "nmSC", "periC"))
Idents(sc_merge_subset) <- sc_merge_subset$cluster

dotPlot(
  path = file.path("lookup", "markers.csv"),
  object = sc_merge_subset,
  par = "novel",
  dot_min = 0.07,
  height = 2,
  width = 4,
)

pns_sn_sciatic_milbrandt_subset <- subset(pns_sn_sciatic_milbrandt, subset = cluster %in% c("mySC", "nmSC", "PnC"))
Idents(pns_sn_sciatic_milbrandt_subset) <- factor(pns_sn_sciatic_milbrandt_subset$cluster, levels = c("mySC", "nmSC", "PnC"))

dotPlot(
  path = file.path("lookup", "markers.csv"),
  object = pns_sn_sciatic_milbrandt_subset,
  par = "novel",
  dot_min = 0.07,
  height = 2,
  width = 4,
  ortho = "human2mouse"
)
