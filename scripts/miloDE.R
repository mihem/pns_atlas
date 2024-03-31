# libraries
library(miloR)
library(miloDE)
library(SingleCellExperiment())
library(Seurat)
library(qs)
library(tidyverse)
library(BiocParallel)
library(BPCells)

# read preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 6)
sc_merge$level0 <- ifelse(sc_merge$level2 == "CTRL", "CTRL", "PNP")

remotes::install_github(("mihem/miloDE@faster_plotting"))
detach(package:miloDE, unload = TRUE)

# # convert to sparse matrix and scale data, necessary for miloDE
# however this needs a lot of memory, so first use on disk matrix
# and then later add on memory matrix to milo DE object
# sc_merge$RNA$counts <- as(object = sc_merge[["RNA"]]$counts, Class = "dgCMatrix")
# sc_merge$RNA$data <- as(object = sc_merge[["RNA"]]$data, Class = "dgCMatrix")

# downsample ---
sc_merge_small <- subset(sc_merge, downsample = 1000)
sc_merge_small <- ScaleData(sc_merge_small)
qsave(sc_merge_small, file.path("objects", "sc_merge_small.qs"))

# convert to sparse matrix and scale data
sc_merge_small$RNA$counts <- as(object = sc_merge_small[["RNA"]]$counts, Class = "dgCMatrix")
sc_merge_small$RNA$data <- as(object = sc_merge_small[["RNA"]]$data, Class = "dgCMatrix")

# remove unnecessary assays and dims
sc_diet <- DietSeurat(
    sc_merge,
    # sc_merge_small,
    counts = TRUE,
    data = TRUE,
    scale.data = TRUE,
    assays = "RNA",
    dimreducs = c("integrated.scvi.full", "umap.scvi.full")
)

#sanity check
identical(sc_diet$RNA$counts, sc_merge_small$RNA$counts)
identical(sc_diet$RNA$data, sc_merge_small$RNA$data)

# very important to set assay when using as.SingleCellExperiment
# use RNA because it holds counts + logcounts (data)
# integrated does not hold counts and only has 2k top variable genes
sce <- as.SingleCellExperiment(sc_diet, assay = "RNA")
counts(sce)
logcounts(sce)

#  mulit core for miloDE
ncores <- 6
mcparam <- MulticoreParam(workers = ncores)
register(mcparam)

# for integrated use PCA (corrected PCA based on integration assay)
# for harmony use harmony
# for scvi use integrated.scvi
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

qs::qsave(milo_DE, file.path("objects", "milo_DE.qs"))
# qs::qsave(milo_DE, file.path("objects", "milo_DE_downsampled.qs"))
milo_DE <- qs::qread(file.path("objects", "milo_DE.qs"))

nhoods_sce <- miloR::nhoods(milo_DE)
nhood_stat_ct <- data.frame(Nhood = 1:ncol(nhoods_sce), Nhood_center = colnames(nhoods_sce))
nhood_stat_ct <- miloR::annotateNhoods(milo_DE, nhood_stat_ct, coldata_col = "cluster")

milo_hood_plot <-
  miloDE::plot_milo_by_single_metric(
    milo_DE,
    nhood_stat_ct,
    colour_by = "cluster",
    layout = "UMAP.SCVI.FULL",
    size_range = c(1.5, 3),
    edge_width = c(0.01, 0.05)
  ) +
  scale_fill_manual(values = my_cols_25, name = "cluster")

ggsave(plot = milo_hood_plot, file.path("results", "miloDE", "milo_rna_nhood.pdf"), width = 10, height = 7)
ggsave(file.path("results", "miloDE", "milo_rna_nhood_downsampled.pdf"), width = 10, height = 7)


## #calculate AUC
## stat_auc <- miloDE::calc_AUC_per_neighbourhood(
##   milo_DE,
##   sample_id = "sample",
##   condition_id = "disease",
##   min_n_cells_per_sample = 1,
##   BPPARAM = mcparam
## )

## qs::qsave(stat_auc, "milo_stat_auc.qs")

## plot_milo_by_single_metric(
##   milo_DE,
##   stat_auc,
##   colour_by = "auc",
##   layout = "UMAP",
##   size_range = c(1.5,3),
##   edge_width = c(0.2,0.5)
## ) +
##   scale_fill_viridis(name = "AUC")
## ggsave(file.path("miloDE", "AUC_dura_small.pdf"), width = 7, height = 7)

milo_DE <- qs::qread(file.path("objects", "milo_DE.qs"), nthread = 6)

counts(miloDE) <- as(object = sc_merge[["RNA"]]$counts, Class = "dgCMatrix")
qsave(milo_DE, file.path("objects", "milo_DE.qs"))

str(counts(miloDE), max.level = 3)
str(logcounts(miloDE), max.level = 3)
scMisc::lss()

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
    contrasts = c("level2CIDP - level2CTRL"),
    # contrasts = c("level0PNP - level0CTRL"),
    # contrasts = c("level2VN - level2CTRL"),
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


# system.time(
#   length(miloR:::graph(milo_DE))
# )

system.time(
p1 <- plot_milo_by_single_metric(
  milo_DE,
  stat_de_magnitude,
  colour_by = "n_DE_genes",
  layout = "UMAP.SCVI.FULL",
  # size_range = c(0.2, 3),
  size_range = c(0.5, 5),
  edge_width = c(0.1, 1.0), 
  edge_weight.thres = 10
) +
  viridis::scale_fill_viridis(name = "# DE genes", option = "inferno")
)

ggsave(plot = p1, filename = file.path("results", "miloDE", "milo_DE_PNP_CTRL_v3.pdf"), width = 6, height = 6, device = cairo_pdf)
ggsave(plot = p1, filename = file.path("results", "miloDE", "milo_DE_VN_CTRL.pdf"), width = 6, height = 6, device = cairo_pdf)
ggsave(plot = p1, filename = file.path("results", "miloDE", "milo_DE_CIDP_CTRL.pdf"), width = 6, height = 6, device = cairo_pdf)
ggsave(plot = p1, filename = file.path("results", "miloDE", "milo_DE_CIAP_CTRL.pdf"), width = 6, height = 6, device = cairo_pdf)

# p2 <- plot_milo_by_single_metric(
#   milo_DE,
#   stat_de_magnitude,
#   colour_by = "n_specific_DE_genes",
#   layout = "UMAP.SCVI.FULL",
#   size_range = c(0.2, 3),
#   ## edge_width = c(0.01, 0.1)
#   edge_weight.thresh = 10
# ) +
#   viridis::scale_fill_viridis(name = "# specific\nDE genes", option = "inferno")

# p_gene <- plot_DE_single_gene(
#   milo_DE,
#   de_stat = de_stat, 
#   gene = "PPIAL4G",
#   layout = "UMAP.SCVI.FULL",
#   set_na_to_0 = FALSE,
#   edge_weight.thresh = 10,
#   alpha = 1
# ) +
#   viridis::scale_fill_viridis(name = "# specific\nDE genes", option = "inferno")

# ggsave(plot = p3, file.path("results", "miloDE", "milo_DE_PNP_CTRL_IFI44L.pdf"), width = 6, height = 6, device = cairo_pdf)
# p1_p2 <- patchwork::wrap_plots(p1, p2)




# stat_de_magnitude |>
#   # dplyr::arrange(desc(n_specific_DE_genes)) |>
#   dplyr::arrange(desc(n_DE_genes)) |>
#   head()

# p_genes_37 <-
#   de_stat@assays@data@listData$pval_corrected_across_genes |>
#   data.frame() |>
#   select(X37) |>
#   rownames_to_column(var = "gene") |>
#   tibble() |>
#   dplyr::filter(X37 < 0.1)

# assay_de_specific <- assay(de_stat, "pval_corrected_across_nhoods")
# idx_de_specific <- which(is.na(assay_de_specific))
# assay_de_specific[idx_de_specific] <- 1
# assay_de_specific_z <- t(apply(assay_de_specific, 1, function(x) {(x - mean(x, na.rm = TRUE))/sd(x, na.rm  = TRUE)}))
# head(sort(assay_de_specific_z[,37][assay_de_specific_z[,37] < -3]), 25)

# assay_de_specific["PPIAL4G",]
# assay_de_specific_z["HLA-DPB1",]


# assay_de_genes <- assay(de_stat, "pval_corrected_across_genes")

# idx_de_genes <- which(is.na(assay_de_genes))
# assay_de_genes[idx_de_genes] <- 1
# assay_de_genes_z <- t(apply(assay_de_genes, 1, function(x) {(x - mean(x, na.rm = TRUE))/sd(x, na.rm  = TRUE)}))

# head(sort(assay_de_genes_z[,425][assay_de_genes_z[,425] < -3]), 25)


# fp_1 <- FeaturePlot(sc_merge, reduction = "umap.scvi.full", features = c("PPIAL4G"), split.by = "level0", pt.size = .01, raster = FALSE, order = TRUE)
# ggsave(plot = fp_1, file.path("results", "miloDE", "milo_DE_PPIAL4G.png"), width = 12, height = 7)

# fp_2 <- FeaturePlot(sc_merge, reduction = "umap.scvi.full", features = c("HLA-DPB1"), split.by = "level0", pt.size = .01, raster = FALSE, order = TRUE)
# ggsave(plot = fp_2, file.path("results", "miloDE", "milo_DE_HLA-DPB1.png"), width = 12, height = 7)

# fp_3 <- FeaturePlot(sc_merge, reduction = "umap.scvi.full", features = c("SERPINA1"), split.by = "level0", pt.size = .01, raster = FALSE, order = TRUE)
# ggsave(plot = fp_3, file.path("results", "miloDE", "milo_SERPINA1.png"), width = 12, height = 7)
