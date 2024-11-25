#===============================================================================
# Single-Cell Data Preparation and Analysis using Milo
#===============================================================================
# Purpose: Prepare single-cell data for differential expression analysis using 
# the Milo package, including data conversion, scaling, and neighborhood assignment.
#
# Methods: 
# - Load preprocessed single-cell data
# - Convert data to sparse matrix format
# - Scale data and prepare for Milo analysis
# - Assign neighborhoods and perform differential expression analysis
#===============================================================================

# Libraries and Setup ----
library(miloR)
library(miloDE)
library(SingleCellExperiment)
library(Seurat)
library(qs)
library(tidyverse)
library(BiocParallel)
library(BPCells)

# Load Preprocessed Data ----
# Load the preprocessed single-cell data
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 6)
sc_merge$level0 <- ifelse(sc_merge$level2 == "CTRL", "CTRL", "PNP")

# Data Conversion and Scaling ----
# Convert to sparse matrix and scale data
sc_merge$RNA$counts <- as(object = sc_merge[["RNA"]]$counts, Class = "dgCMatrix")
sc_merge$RNA$data <- as(object = sc_merge[["RNA"]]$data, Class = "dgCMatrix")

# Remove unnecessary assays and dimensions
sc_diet <- DietSeurat(
    sc_merge,
    counts = TRUE,
    data = TRUE,
    scale.data = TRUE,
    assays = "RNA",
    dimreducs = c("integrated.scvi.full", "umap.scvi.full")
)

# Sanity check for data integrity
identical(sc_diet$RNA$counts, sc_merge$RNA$counts)
identical(sc_diet$RNA$data, sc_merge$RNA$data)

# Convert to SingleCellExperiment
sce <- as.SingleCellExperiment(sc_diet, assay = "RNA")

# Neighborhood Assignment ----
# Assign neighborhoods using Milo
set.seed(123)
milo_DE <- assign_neighbourhoods(
    sce,
    k = 30,
    prop = 0.1,
    d = 30,
    order = 2,
    filtering = TRUE,
    reducedDim_name = "INTEGRATED.SCVI.FULL",
    verbose = TRUE
)

# Save Milo object
qs::qsave(milo_DE, file.path("objects", "milo_DE.qs"))

# Neighborhood Annotation and Plotting ----
# Annotate neighborhoods and create plots
nhoods_sce <- miloR::nhoods(milo_DE)
nhood_stat_ct <- data.frame(Nhood = 1:ncol(nhoods_sce), Nhood_center = colnames(nhoods_sce))
nhood_stat_ct <- miloR::annotateNhoods(milo_DE, nhood_stat_ct, coldata_col = "cluster")

milo_hood_plot <- miloDE::plot_milo_by_single_metric(
    milo_DE,
    nhood_stat_ct,
    colour_by = "cluster",
    layout = "UMAP.SCVI.FULL",
    size_range = c(1.5, 3),
    edge_width = c(0.01, 0.05)
) + scale_fill_manual(values = my_cols_25, name = "cluster")

ggsave(plot = milo_hood_plot, file.path("results", "miloDE", "milo_rna_nhood.pdf"), width = 10, height = 7)

# Differential Expression Analysis ----
# Perform differential expression analysis using Milo
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

# Save Results ----
# Save the differential expression statistics
qs::qsave(de_stat, file.path("objects", "milo_de_stat_vn_ctrl.qs"))
qs::qsave(de_stat, file.path("objects", "milo_de_stat_pnp_ctrl.qs"))
qs::qsave(de_stat, file.path("objects", "milo_de_stat_vn_ctrl_downsampled.qs"))
qs::qsave(de_stat, file.path("objects", "milo_de_stat_pnp_ctrl_downsampled.qs"))


## system("systemctl suspend")
milo_DE <- qs::qread(file.path("objects", "milo_DE.qs"), nthread = 4)
de_stat <- qs::qread(file.path("objects", "milo_de_stat_pnp_ctrl.qs"), nthread = 6)
de_stat <- qs::qread(file.path("objects", "milo_de_stat_vn_ctrl.qs"), nthread = 6)
de_stat <- qs::qread(file.path("objects", "milo_de_stat_cidp_ctrl.qs"), nthread = 6)
de_stat <- qs::qread(file.path("objects", "milo_de_stat_ciap_ctrl.qs"), nthread = 6)

stat_de_magnitude <- rank_neighbourhoods_by_DE_magnitude(de_stat)

system.time(
p1 <- plot_milo_by_single_metric(
  milo_DE,
  stat_de_magnitude,
  colour_by = "n_DE_genes",
  layout = "UMAP.SCVI.FULL",
  size_range = c(0.5, 5),
  edge_width = c(0.1, 1.0), 
  edge_weight.thres = 10
) +
  viridis::scale_fill_viridis(name = "# DE genes", option = "inferno")
)

ggsave(plot = p1, filename = file.path("results", "miloDE", "milo_DE_PNP_CTRL.pdf"), width = 6, height = 6, device = cairo_pdf)
ggsave(plot = p1, filename = file.path("results", "miloDE", "milo_DE_VN_CTRL.pdf"), width = 6, height = 6, device = cairo_pdf)
ggsave(plot = p1, filename = file.path("results", "miloDE", "milo_DE_CIDP_CTRL.pdf"), width = 6, height = 6, device = cairo_pdf)
ggsave(plot = p1, filename = file.path("results", "miloDE", "milo_DE_CIAP_CTRL.pdf"), width = 6, height = 6, device = cairo_pdf)
