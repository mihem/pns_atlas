# ===============================================================================
# Figure Reproducibility Preparation Script
# ===============================================================================
# Purpose: Prepare Seurat object for qmd file (reproducibility the figures)
# ===============================================================================

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

# Remove unnecessary data to further reduce the object size
ic_figure <- DietSeurat(
    ic,
    counts = TRUE,
    data = FALSE,
    scale.data = FALSE,
    assays = "RNA",
    dimreducs = c("umap.rpca")
)

ic_figure$RNA$counts <- NULL
ic_figure$RNA$scale.data <- NULL

# Remove unnecessary data to further reduce the object size
ic_figure$RNA$counts <- NULL
ic_figure$RNA$scale.data <- NULL

ic_figure@meta.data <- data.frame()
ic_figure@commands <- list()
ic_figure@tools <- list()

qsave(ic_figure, file.path("docs", "ic_figure.qs"))
ic_figure <- qread(file.path("docs", "ic_figure.qs"))

# subset the object to only include the B and Plasma clusters
b_plasma <- subset(ic, idents = c("Plasma", "B"))

# Remove unnecessary data to further reduce the object size
b_plasma_figure <- DietSeurat(
    b_plasma,
    counts = TRUE,
    data = FALSE,
    scale.data = FALSE,
    assays = "RNA",
    dimreducs = NULL
)

b_plasma_figure@meta.data <-
    b_plasma_figure@meta.data |>
    dplyr::select(ic_cluster)

b_plasma_figure$RNA$counts <- NULL
b_plasma_figure$RNA$scale.data <- NULL
b_plasma_figure@commands <- list()
b_plasma_figure@tools <- list()

qsave(b_plasma_figure, file.path("docs", "b_plasma_figure.qs"))

# save enrichment results as qs object
enrichr_macro18 <- readxl::read_excel(file.path("results", "enrichr", "enrichr_Macro18.xlsx"), sheet = "GO_Biological_Process_2023")
qsave(enrichr_macro18, file.path("docs", "enrichr_macro18.qs"))

# Figure 3 ----
propeller_PNP_CTRL <-
    scMisc::propellerCalc(
        seu_obj1 = sc_merge,
        condition1 = "PNP",
        condition2 = "CTRL",
        cluster_col = "cluster",
        meta_col = "level0",
        lookup = sample_lookup,
        sample_col = "sample",
        formula = "~0 + level0",
        min_cells = 30
    ) |>
    dplyr::filter(abs(log2ratio) > 0.5)

qsave(propeller_PNP_CTRL, file.path("docs", "propeller_PNP_CTRL.qs"))

propeller_PNP_CTRL_ic <-
    scMisc::propellerCalc(
        seu_obj1 = ic,
        condition1 = "PNP",
        condition2 = "CTRL",
        cluster_col = "ic_cluster",
        meta_col = "level0",
        lookup = sample_lookup,
        sample_col = "sample",
        formula = "~0 + level0",
        min_cells = 30
    ) |>
    dplyr::filter(abs(log2ratio) > 0.5)

qsave(propeller_PNP_CTRL_ic, file.path("docs", "propeller_PNP_CTRL_ic.qs"))

# Figure 3B CNA  ----
# PNP
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)
sc_merge$pnp <- as.numeric(sc_merge$level2 != "CTRL")
sc_merge$sex_numeric <- as.numeric(sc_merge$sex == "male")
sc_merge$center_numeric <- as.numeric(factor(sc_merge$center))

cna_pnp <- association.Seurat(
    seurat_object = sc_merge, 
    test_var = 'pnp', 
    samplem_key = 'sample', 
    graph_use = 'RNA_nn', 
    verbose = TRUE,
    batches = NULL, ## no batch variables to include, only works with matched design https://github.com/immunogenomics/cna/issues/11
    covs = c("sex_numeric", "age")
)
cna_pnp$cna_ncorrs_pnp <- cna_pnp$cna_ncorrs

# g ratio
cna_gratio <- association.Seurat(
    seurat_object = sc_merge, 
    test_var = "g_ratio", 
    samplem_key = 'sample', 
    graph_use = 'RNA_nn', 
    verbose = TRUE,
    batches = NULL, ## no batch variables to include
    covs = c("sex_numeric", "age") 
)

cna_pnp$cna_ncorrs_gratio <- cna_gratio$cna_ncorrs

# Remove unnecessary data to further reduce the object size
cna_pnp_gratio_figure <- DietSeurat(
    cna_pnp,
    counts = TRUE,
    data = FALSE,
    scale.data = FALSE,
    assays = "RNA",
    dimreducs = c("umap.scvi.full")
)

cna_pnp_gratio_figure$RNA$counts <- NULL
cna_pnp_gratio_figure$RNA$scale.data <- NULL
cna_pnp_gratio_figure@commands <- list()
cna_pnp_gratio_figure@tools <- list()
cna_pnp_gratio_figure@meta.data <-
    cna_pnp_gratio_figure@meta.data |>
    tibble::rownames_to_column("barcode") |>
    dplyr::select(barcode, cna_ncorrs_pnp, cna_ncorrs_gratio) |>
    tibble::column_to_rownames("barcode")

qsave(cna_pnp_gratio_figure, file.path("docs", "cna_pnp_gratio_figure.qs"))

