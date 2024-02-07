# libraries ----
library(miloR)
library(SingleCellExperiment())
library(Seurat)
library(qs)
library(tidyverse)
library(BiocParallel)

# read preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)
sc_merge <- ScaleData(sc_merge)

# no necessary for miloR, works with BPCellsMatrix (on disk)
# # convert to sparse matrix and scale data
# sc_merge$RNA$counts <- as(object = sc_merge[["RNA"]]$counts, Class = "dgCMatrix")
# sc_merge$RNA$data <- as(object = sc_merge[["RNA"]]$data, Class = "dgCMatrix")
# sc_merge <- ScaleData(sc_merge)
# sc_merge <- ScaleData(sc_merge)

scMisc::lss()

# downsample ---
sc_merge_small <- subset(sc_merge, downsample = 1000)
sc_merge_small <- ScaleData(sc_merge_small)
qsave(sc_merge_small, file.path("objects", "sc_merge_small.qs"))

# convert to sparse matrix and scale data
sc_merge_small$RNA$counts <- as(object = sc_merge_small[["RNA"]]$counts, Class = "dgCMatrix")
sc_merge_small$RNA$data <- as(object = sc_merge_small[["RNA"]]$data, Class = "dgCMatrix")
sc_merge_small <- ScaleData(sc_merge_small)

str(sc_merge_small$RNA$counts)
str(sc_merge_small$RNA$data)
str(sc_merge_small$RNA$scale.data)

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
identical(sc_diet$RNA$counts, sc_merge$RNA$counts)
identical(sc_diet$RNA$data, sc_merge$RNA$data)

# somehow counts data are not kept in the downsampled object
sc_diet$RNA$data <- sc_merge_small$RNA$data

# very important to set assay when using as.SingleCellExperiment
# use RNA because it holds counts + logcounts (data)
# integrated does not hold counts and only has 2k top variable genes
# for integrated use PCA (corrected PCA based on integration assay)
# for harmony use harmony
# for scvi use integrated.scvi
sce <- as.SingleCellExperiment(sc_diet, assay = "RNA")

scMisc::lss()


#  mulit core for miloDE
ncores <- 6
mcparam <- MulticoreParam(workers = ncores)
register(mcparam)

milo_obj <- Milo(sce)
str(counts(milo_obj)@matrix, max.level = 2)
str(logcounts(milo_obj)@matrix, max.level = 2)

# takes very long with all cells
set.seed(123)
milo_obj <- buildGraph(milo_obj, k = 30, d = 30, reduced.dim = "INTEGRATED.SCVI.FULL")

#very important to use HARMONY for reduced_dims, if Seurat CCA is run use PCA, if SCVI use e.g. INTEGTRATED:SCVI.FULL
set.seed(123)

milo_obj <- makeNhoods(
  milo_obj,
#   prop = 0.1,
   prop = 0.05,
  k = 30,
  ## k = 50,
  d = 30,
  refined = TRUE,
  reduced_dims = "INTEGRATED.SCVI.FULL",
  refinement_scheme = "graph"
)

plotNhoodSizeHist(milo_obj)

colData(milo_obj)

milo_obj <- miloR::countCells(milo_obj, meta.data = data.frame(colData(milo_obj)), sample = "sample")

head(miloR::nhoodCounts(milo_obj))

# defining experimental design ------
sc_design <- data.frame(colData(milo_obj))[,c("sample", "level2", "sex", "age", "center")]
sc_design <- dplyr::distinct(sc_design)
rownames(sc_design) <- sc_design$sample

# test DA ----
# design changes the p values drastically
da_results <- miloR::testNhoods(
    milo_obj,
    # design = ~ 0 + level2 + sex + age + center,
    design = ~ 0 + level2,
    design.df = sc_design,
    model.contrasts = c("level2VN - level2CTRL"),
    ## model.contrasts = c("aie_typeCASPR2 - aie_typeIIH"),
    ## model.contrasts = c("aie_typeCASPR2 - aie_typeLGI1"),
    fdr.weighting = "graph-overlap",
    ## fdr.weighting = "k-distance",
    reduced.dim = "INTEGRATED.SCVI.FULL"
)

# inpsect DA results -----
da_results |>
  arrange(SpatialFDR) |>
  head()

ggplot(da_results, aes(logFC, -log10(SpatialFDR))) +
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)

milo_obj <- miloR::buildNhoodGraph(milo_obj)

nh_graph_pl <- miloR::plotNhoodGraphDA(
  milo_obj,
  da_results,
  layout = "UMAP.SCVI.FULL",
  alpha = .1,
  size_range = c(0.2, 2),
  node_stroke = 0.1,
)

ggsave(plot = nh_graph_pl, file.path("results", "miloR", "neighborhood_VN_CTRL.pdf"), width = 7, height = 7)
ggsave(plot = nh_graph_pl, file.path("results", "miloR", "neighborhood_VN_CTRL_all.png"), width = 10, height = 10)

da_results <- miloR::annotateNhoods(milo_obj, da_results, coldata_col = "cluster")

## ggplot(da_results, aes(cluster_fraction)) + geom_histogram(bins = 50)
str(da_results)

da_results$cluster_safe <- ifelse(da_results$cluster_fraction < 0.7, "Mixed", da_results$cluster)

miloR::plotDAbeeswarm(da_results, group.by = "cluster_safe", alpha = 0.1)

ggsave(file.path("results", "miloR", "beeswarm_VN_CTRL_all.pdf"), width = 7, height = 7)

qs::qsave(milo_obj, file.path("objects", "milo_object_VN_CTRL_all.qs"))
qs::qsave(da_results, file.path("objects", "milo_da_results_VN_CTRL_all.qs"))
