#===============================================================================
# Documented: Prepare Seurat object for mrvi analysis and plot mrvi results
# after mrvi.py has been run
#===============================================================================
# Libraries ----
# Load necessary libraries for data processing and visualization
library(Seurat)
library(BPCells)
library(SeuratObject)
library(tidyverse)
library(qs)
library(sceasy)
library(pheatmap)
library(pals)
library(viridis)

# Read preprocessed data ----
# Load the merged single-cell dataset
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)

# Downsample data ----
# Reduce dataset size for faster processing
sc_merge_small <- subset(sc_merge, downsample = 1000)

# Remove unnecessary assays and dimensions ----
# Simplify the Seurat object by keeping only relevant assays and dimensional reductions
sc_diet <- DietSeurat(
    sc_merge_small,
    counts = TRUE,
    data = TRUE,
    scale.data = TRUE,
    assays = "RNA",
    dimreducs = c("integrated.scvi.full", "umap.scvi.full")
)

# Prepare data for Python ----
# Convert cluster identifiers to character for compatibility
sc_diet$cluster <- as.character(Idents(sc_diet))
DimPlot(sc_diet, label = TRUE, group.by = "cluster")

# Convert Seurat v5 to Seurat v3 format ----
sc_diet[["RNA3"]] <- as(object = sc_diet[["RNA"]], Class = "Assay")
DefaultAssay(sc_diet) <- "RNA3"
sc_diet[["RNA"]] <- NULL
sc_diet <- RenameAssays(object = sc_diet, RNA3 = 'RNA')

# Convert to AnnData format using sceasy ----
sceasy::convertFormat(sc_diet, from = "seurat", to = "anndata", outFile = file.path("objects", "sc_diet.h5ad"))

################################################################################################################
# after mrvi.py has binished contiue with this script

# MRVI Analysis ----
# Load MRVI results for all clusters combined
mrvi_all_average <-
    read.csv(file.path("results", "mrvi", "mrvi_average_all.csv")) |>
    tibble::column_to_rownames("X")

# Prepare metadata ----
# Merge sample metadata with histological measurements
histo_lookup <-
    tibble(
        sample = sc_merge$sample,
        axon_normal = sc_merge$axon_normal,
        log_axon_normal = sc_merge$log_axon_normal,
        g_ratio = sc_merge$g_ratio,
        axon_diameter = sc_merge$axon_diameter
    ) |>
    distinct() 

# only for CIDP, CIAP, CTRL and VN ----
metadata <- 
    read_csv(file.path("lookup", "sample_lookup.csv")) |>
    left_join(histo_lookup, join_by(sample == sample)) |>
    dplyr::filter(level2 %in% c("CIDP", "CIAP", "CTRL", "VN")) |>
    select(sample, level2, INCAT, log_axon_normal, axon_diameter, g_ratio) |>
    rename(log_axon_count = log_axon_normal) |>
    mutate(level0 = if_else(level2 == "CTRL", "CTRL", "PNP")) |>
    relocate(level0, level2, INCAT, log_axon_count, axon_diameter, g_ratio) |>
    data.frame()

mrvi_data <- 
    mrvi_all_average |>
    rownames_to_column("sample") |>
    left_join(metadata, join_by(sample == sample)) |>
    dplyr::filter(level2 %in% c("CIDP", "CIAP", "CTRL", "VN")) |>
    column_to_rownames("sample") |>
    select(metadata$sample)

# make sure that rownames and columnnames of data and metadata are the same
rownames(metadata) <- colnames(mrvi_data)

# create colors for annotations
phmap_cols_level2 <- sc_merge@misc$level2_cols[c("CIAP", "CIDP", "CTRL", "VN")]

# plot colors of set1
# phmap_cols_level0 <- RColorBrewer::brewer.pal(8, "Set1")[c(7, 6)]
phmap_cols_level0 <- pals::cols25()[c(7, 6)]
names(phmap_cols_level0) <- unique(metadata[["level0"]])

phmap_cols_INCAT <- c("white", RColorBrewer::brewer.pal(length(unique(metadata$INCAT))-1, "Reds"))
names(phmap_cols_INCAT) <- c("-", "1", "2", "3", "4", "5", "6")

phmap_cols_axon_count <- viridis::magma(length(unique(metadata$log_axon_count)))
names(phmap_cols_axon_count) <- unique(metadata[["log_axon_count"]])

phmap_cols_axon_diameter <- viridis::magma(length(unique(metadata$axon_diameter)))
names(phmap_cols_axon_diameter) <- unique(metadata$axon_diameter)

phmap_cols_gratio <- viridis::magma(length(unique(metadata$g_ratio)))
names(phmap_cols_gratio) <- unique(metadata$g_ratio)

phmap_list <-
    list(
        level0 = phmap_cols_level0,
        level2 = phmap_cols_level2,
        INCAT = phmap_cols_INCAT,
        log_axon_count = phmap_cols_axon_count,
        axon_diameter = phmap_cols_axon_diameter,
        g_ratio = phmap_cols_gratio
    )

# Generate heatmaps for MRVI data ----
phmap_mrvi <-
    pheatmap(
        mrvi_data,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        color = viridis::magma(100),
        cellwidth = 10,
        cellheight = 10,
        treeheight_row = 15,
        treeheight_col = 15,
        clustering_distance_rows = "euclidean",
        clustering_distance_cols = "euclidean",
        clustering_method = "ward.D2",
        border_color = NA,
        cutree_cols = 5,
        cutree_rows = 5,
        annotation_row = select(metadata, !sample),
        annotation_colors = phmap_list,
        main = "All"
    )

# Save heatmap as PDF ----
pdf(file.path("results", "mrvi", "heatmap_mrvi_average_all.pdf"), width = 10, height = 10)
print(phmap_mrvi)
dev.off()
