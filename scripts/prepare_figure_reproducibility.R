# ===============================================================================
# Figure Reproducibility Preparation Script
# ===============================================================================
# Purpose: Prepare Seurat object for qmd file (reproducibility the figures)
# ===============================================================================

# Load necessary libraries
library(Seurat)
library(qs)
library(tidyverse)
library(miloDE)
library(SingleCellExperiment)

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

# Prepare cna_incat
sc_merge_pnp <- subset(sc_merge, subset = level2 %in% c("CTRL"), invert = TRUE)
sc_merge_pnp$incat_numeric <- as.numeric(sc_merge_pnp$incat)
sc_merge$sex_numeric <- as.numeric(sc_merge$sex == "male")

cna_incat <- association.Seurat(
    seurat_object = sc_merge_pnp, 
    test_var = 'incat_numeric', 
    samplem_key = 'sample', 
    graph_use = 'RNA_nn', 
    verbose = TRUE,
    batches = NULL, 
    covs = c("sex_numeric", "age")
)

# Remove unnecessary data to further reduce the object size
cna_incat_figure <- DietSeurat(
    cna_incat,
    counts = TRUE,
    data = FALSE,
    scale.data = FALSE,
    assays = "RNA",
    dimreducs = c("umap.scvi.full")
)

cna_incat_figure$RNA$counts <- NULL
cna_incat_figure$RNA$scale.data <- NULL
cna_incat_figure@commands <- list()
cna_incat_figure@tools <- list()
cna_incat_figure@meta.data <-
    cna_incat_figure@meta.data |>
    tibble::rownames_to_column("barcode") |>
    dplyr::select(barcode, cna_ncorrs) |>
    tibble::column_to_rownames("barcode")

qsave(cna_incat_figure, file.path("docs", "cna_incat_figure.qs"))

# Figure 3C ----
sheets <- readxl::excel_sheets(path = file.path("results", "de", "pnp_ctrl_pseudobulk.xlsx"))
cl_sig <-
        lapply(
            sheets,
            function(sheet) {
                readxl::read_xlsx(path = file.path("results", "de", "pnp_ctrl_pseudobulk.xlsx"), sheet = sheet) |>
                    dplyr::filter(p_val_adj < 0.05) |>
                    nrow()
            }
        )
result <- tibble(
        cluster = sheets,
        n = unlist(cl_sig),
    )

qs::qsave(result, file.path("docs", "pnp_ctrl_pseudobulk_de.qs"))

# Figure 3D ----
milo_DE <- qs::qread(file.path("objects", "milo_DE.qs"), nthread = 4)

# Print size of each slot in the Milo object
slots <- slotNames(milo_DE)
for(slot in slots) {
    size <- object.size(slot(milo_DE, slot))
    print(sprintf("Slot %s: %s", slot, format(size, units = "auto")))
}

de_stat <- qs::qread(file.path("objects", "milo_de_stat_pnp_ctrl.qs"), nthread = 6)
stat_de_magnitude <- rank_neighbourhoods_by_DE_magnitude(de_stat)

milo_DE <- milo_DE_backup
milo_DE_backup <- milo_DE

# Remove unnecessary data to reduce the object size
empty_matrix <- Matrix::Matrix(0, 
                             nrow = nrow(milo_DE), 
                             ncol = ncol(milo_DE), 
                             sparse = TRUE)
rownames(empty_matrix) <- rownames(milo_DE)
colnames(empty_matrix) <- colnames(milo_DE)

for(i in assayNames(milo_DE)) {
    assay(milo_DE, i) <- empty_matrix
}

colData(milo_DE) <- NULL
milo_DE@graph <- list()
reducedDims(milo_DE) <- reducedDims(milo_DE)["UMAP.SCVI.FULL"]

# Print the size of the object
print(object.size(milo_DE), units = "Gb")

milo_figure <- list(
    obj = milo_DE,
    stat = stat_de_magnitude
)

qsave(milo_figure, file.path("docs", "milo_figure.qs"))

# Figure 3E ----
# calculate DE PNP vs Ctrl for each cluster
dePseudo <- function(seu_obj, cell_type_col, label_col) {
    res <- Libra::run_de(seu_obj,
        replicate_col = "sample",
        cell_type_col = cell_type_col,
        label_col = label_col,
        min_cells = 3,
        min_reps = 2,
        min_feature = 0,
        de_family = "pseudobulk",
        de_method = "edgeR",
        de_type = "LRT",
        n_threads = 6
    )
    res <- arrange(res, cell_type, desc(avg_logFC))
    res_split <- split(res, res$cell_type)
    return(res_split)
}

de_pseudo_pnp_ctrl <- dePseudo(sc_merge, cell_type_col = "cluster", label_col = "level0")
qsave(de_pseudo_pnp_ctrl, file.path("docs", "de_pseudo_pnp_ctrl.qs"))

qsave(de_pseudo_pnp_ctrl, file.path("docs", "de_pseudo_pnp_ctrl.qs"))
