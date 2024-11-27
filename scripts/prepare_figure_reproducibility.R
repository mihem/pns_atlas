#===============================================================================
# Figure Reproducibility Preparation Script
#===============================================================================
# Purpose: Prepare Seurat object for the qmd file (reproducibility the figures)
#===============================================================================

# Load necessary libraries
library(Seurat)
library(qs)

# Figure 1 ----
# Prepare Seurat object for reproducibility
# DietSeurat reduces the size of the Seurat object by keeping only the necessary data
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)
umap_figure <- DietSeurat(
    sc_merge,
    counts = TRUE,
    data = FALSE,
    scale.data = FALSE,
    assays = "RNA",
    dimreducs = c("umap.scvi.full")
)

# Remove unnecessary data to further reduce the object size
umap_figure$RNA$counts <- NULL
umap_figure$RNA$data <- NULL
umap_figure$RNA$scale.data <- Matrix::Matrix(
  0,
  nrow = nrow(umap_figure$RNA),
  ncol = ncol(umap_figure$RNA)
)

umap_figure@meta.data <- data.frame()
umap_figure@commands <- list()
umap_figure@tools <- list()

# Save the processed Seurat object to a file
qs::qsave(umap_figure, file.path("docs", "umap_figure.qs"))

# Figure 2 ----
ic <- qs::qread(file.path("objects", "ic.qs"), nthread = 4)

ic_umap_figure <- DietSeurat(
    ic,
    counts = TRUE,
    data = FALSE,
    scale.data = FALSE,
    assays = "RNA",
    dimreducs = c("umap.rpca")
)

# Remove unnecessary data to further reduce the object size
ic_umap_figure$RNA$counts <- NULL
ic_umap_figure$RNA$data <- NULL
ic_umap_figure$RNA$scale.data <- Matrix::Matrix(
  0,
  nrow = nrow(ic_umap_figure$RNA),
  ncol = ncol(ic_umap_figure$RNA)
)

ic_umap_figure@meta.data <- data.frame()
ic_umap_figure@commands <- list()
ic_umap_figure@tools <- list()

qsave(ic_umap_figure, file.path("docs", "ic_umap_figure.qs"))

scMisc::lss()
str(ic_umap_figure, max.level = 2)
print(object.size(ic_umap_figure@tools), units = "Mb")


