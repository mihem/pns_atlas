# libraries
library(miloR)
library(miloDE)
library(SingleCellExperiment())
library(qs)
library(tidyverse)
library(BiocParallel)
library(viridis)
library(Seurat)
library(BPCells)

print("Loaded libraries")

my_cols_25 <- pals::cols25()

# read preprocessed data ----
ncores <- 6
mcparam <- MulticoreParam(workers = ncores)
register(mcparam)

# read preprocessed data ----
milo_DE <- qs::qread("milo_DE.qs", nthread = 4)
print("Loaded milo_DE object")

# and ~72h for ~400k cells with 1 core)
system.time(
  de_stat <- miloDE::de_test_neighbourhoods(
    milo_DE,
    sample_id = "sample",
    ## design = ~0 + level2 + sex + age + center,
    design = ~ 0 + level2,
    # design = ~ 0 + level0,
    covariates = c("level2"),
    # covariates = c("level0"),
    # contrasts = c("level2VN - level2CTRL"),
    # contrasts = c("level0PNP - level0CTRL"),
    contrasts = c("level2CIDP - level2CTRL"),
    ## subset_nhoods = stat_auc$Nhood[!is.na(stat_auc$auc)],
    output_type = "SCE",
    plot_summary_stat = TRUE,
    layout = "UMAP.SCVI.FULL",
    BPPARAM = NULL,
    # BPPARAM = mcparam,
    verbose = TRUE,
    min_count = 10
  )
)

print("Done with de_test_neighbourhoods")

qs::qsave(de_stat, file.path("milo_de_stat_vn_ctrl.qs"))