
#===============================================================================
# GEO submission
#
# Purpose: Export metadata for GEO submission
#===============================================================================

# Load required libraries ----
library(Seurat)
library(qs)
library(tibble)

# Data preparation ----
# Load preprocessed data and calculate age for each sample
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)
ic <- qs::qread(file.path("objects", "ic.qs"), nthread = 4)

geo_meta_all <- 
    sc_merge@meta.data |>
    tibble::rownames_to_column("barcode") |>
    dplyr::select(
      barcode,
      sample,
      cluster,
      level2,
      nCount_RNA:incat,
      g_ratio,
      axon_diameter,
      milbrandt_sciatic_label_full,
      milbrandt_sciatic_label_full.score,
      suter_p60_label_full,
      suter_p60_label_full.score,
      rosmap_label,
      rosmap_score
    )

write.csv(geo_meta_all, file = file.path("geo", "metadata_all.csv"), row.names = FALSE)

# immune cells
geo_meta_ic <- 
    ic@meta.data |>
    tibble::rownames_to_column("barcode") |>
    dplyr::select(
      barcode,
      sample,
      ic_cluster,
      level2,
      nCount_RNA:incat,
      g_ratio,
      axon_diameter,
      milbrandt_sciatic_label_full,
      milbrandt_sciatic_label_full.score,
      suter_p60_label_full,
      suter_p60_label_full.score
    )

write.csv(geo_meta_ic, file = file.path("geo", "metadata_ic.csv"), row.names = FALSE)
