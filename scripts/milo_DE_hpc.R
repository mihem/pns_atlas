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

# careful! computationally very expensive (22h runtime with 1 core and 170k cells
# and 38h for ~400k cells with 1 core)
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
    contrasts = c("level2VN - level2CTRL"),
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

stat_de_magnitude <- rank_neighbourhoods_by_DE_magnitude(de_stat)

p1 <- plot_milo_by_single_metric(
  milo_DE,
  stat_de_magnitude,
  colour_by = "n_DE_genes",
  layout = "UMAP.SCVI.FULL",
  size_range = c(0.2, 3),
  ## edge_width = c(0.5, 1.0), # does not work
  edge_weight.thres = 10
) +
  viridis::scale_fill_viridis(name = "# DE genes")
  ## ggraph::scale_edge_width(range = c(0.001, 0.01))


p2 <- plot_milo_by_single_metric(
  milo_DE,
  stat_de_magnitude,
  colour_by = "n_specific_DE_genes",
  layout = "UMAP.SCVI.FULL",
  size_range = c(0.2, 3),
  ## edge_width = c(0.01, 0.1)
  edge_weight.thresh = 10
) +
  viridis::scale_fill_viridis(name = "# specific\nDE genes", option = "inferno")

patchwork::wrap_plots(p1, p2)
ggsave(file.path("results", "miloDE", "milo_DE_VN_CTRL.pdf"), width = 12, height = 6, device = cairo_pdf)