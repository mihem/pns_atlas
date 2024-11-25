#===============================================================================
# Differential Expression Analysis using Milo
#===============================================================================
# Purpose: Perform differential expression analysis on single-cell data using 
# the Milo package, leveraging high-performance computing (HPC) resources.
#
# Methods: 
# - Load preprocessed single-cell data
# - Perform differential expression analysis using Milo
# - Save results for further analysis
#===============================================================================

# Libraries and Setup ----
library(miloR)
library(miloDE)
library(SingleCellExperiment)
library(qs)
library(tidyverse)
library(BiocParallel)
library(viridis)
library(Seurat)
library(BPCells)

print("Loaded libraries")

my_cols_25 <- pals::cols25()

# General Settings ----
# Define number of cores for parallel processing
ncores <- 6
mcparam <- MulticoreParam(workers = ncores)
register(mcparam)

# Load the preprocessed Milo object ---
milo_DE <- qs::qread("milo_DE.qs", nthread = 4)
print("Loaded milo_DE object")

# Differential Expression Analysis ----
# Perform differential expression analysis using Milo
# This step may take a significant amount of time depending on the dataset size
system.time(
  de_stat <- miloDE::de_test_neighbourhoods(
    milo_DE,
    sample_id = "sample",
    design = ~ 0 + level2,
    covariates = c("level2"),
    contrasts = c("level2CIDP - level2CTRL"),
    output_type = "SCE",
    plot_summary_stat = TRUE,
    layout = "UMAP.SCVI.FULL",
    BPPARAM = NULL,
    verbose = TRUE,
    min_count = 10
  )
)

print("Done with de_test_neighbourhoods")

# Save Results ----
# Save the differential expression statistics
qs::qsave(de_stat, file.path("milo_de_stat_vn_ctrl.qs"))