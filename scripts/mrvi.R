# libraries  ----
library(Seurat)
library(BPCells)
library(SeuratObject)
library(tidyverse)
library(qs)
library(sceasy)
library(pheatmap)
library(pals)

# read preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)

# downsample ---
sc_merge_small <- subset(sc_merge, downsample = 1000)

# remove unnecessary assays and dims
sc_diet <- DietSeurat(
    sc_merge_small,
    counts = TRUE,
    data = TRUE,
    scale.data = TRUE,
    assays = "RNA",
    dimreducs = c("integrated.scvi.full", "umap.scvi.full")
)

#python has problems with factors
sc_diet$cluster <- as.character(Idents(sc_diet))
DimPlot(sc_diet, label = TRUE, group.by = "cluster")

# convert seuratv5 to seuratv3 data
sc_diet[["RNA3"]] <- as(object = sc_diet[["RNA"]], Class = "Assay")
DefaultAssay(sc_diet) <- "RNA3"
sc_diet[["RNA"]] <- NULL
sc_diet <- RenameAssays(object = sc_diet, RNA3 = 'RNA')

# convert using sceasy ----
sceasy::convertFormat(sc_diet, from = "seurat", to = "anndata", outFile = file.path("objects", "sc_diet.h5ad"))

# reticulate::py_config()

# mrvvi anaylsis for all cluster combined ---
mrvi_all_average <-
    read.csv(file.path("results", "mrvi", "mrvi_average_all.csv")) |>
    tibble::column_to_rownames("X")

# get metadata
axon_count_lookup <-
    tibble(sample = sc_merge$sample, axon_normal = sc_merge$axon_normal, log_axon_normal = sc_merge$log_axon_normal) |>
    distinct() |>
    mutate()

metadata <- 
    read_csv(file.path("lookup", "sample_lookup.csv")) |>
    left_join(axon_count_lookup) |>
    select(level2, INCAT, center, log_axon_normal) |>
    mutate(level0 = if_else(level2 == "CTRL", "CTRL", "PNP")) |>
    # select(!log_axon_normal) |>
    data.frame()


# make sure that rownames and columnnames of data and metadata are the same
rownames(metadata) <- colnames(mrvi_all_average)

# create colors for annotations
phmap_cols_level2 <- pals::cols25(length(unique(metadata$level2)))
names(phmap_cols_level2) <- unique(metadata[["level2"]])

phmap_cols_level0 <- RColorBrewer::brewer.pal(length(unique(metadata$level0)), "Set1")[1:2]
names(phmap_cols_level0) <- unique(metadata[["level0"]])

phmap_cols_INCAT <- c("white", RColorBrewer::brewer.pal(length(unique(metadata$INCAT))-1, "Reds"))
names(phmap_cols_INCAT) <- c("-", "1", "2", "3", "4", "5", "6")

phmap_cols_center <- RColorBrewer::brewer.pal(length(unique(metadata$center)), "Set2")
names(phmap_cols_center) <- unique(metadata[["center"]])

phmap_cols_axon_normal <- viridis::viridis(length(unique(metadata$log_axon_normal)))
names(phmap_cols_axon_normal) <- unique(metadata[["log_axon_normal"]])

phmap_list <-
    list(
        level0 = phmap_cols_level0,
        level2 = phmap_cols_level2,
        INCAT = phmap_cols_INCAT,
        center = phmap_cols_center,
        log_axon_normal = phmap_cols_axon_normal
    )

phmap_mrvi <-
    pheatmap(
        mrvi_all_average,
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
        annotation_row = metadata,
        annotation_colors = phmap_list,
        main = "All"
    )

pdf(file.path("results", "mrvi", "heatmap_mrvi_average_all.pdf"), width = 10, height = 10)
print(phmap_mrvi)
dev.off()

# only for CIDP, CIAP, CTRL and VN ----
metadata <- 
    read_csv(file.path("lookup", "sample_lookup.csv")) |>
    left_join(axon_count_lookup, join_by(sample == sample)) |>
    dplyr::filter(level2 %in% c("CIDP", "CIAP", "CTRL", "VN")) |>
    select(sample, level2, INCAT, center, log_axon_normal) |>
    mutate(level0 = if_else(level2 == "CTRL", "CTRL", "PNP")) |>
    relocate(log_axon_normal, INCAT, center, level2, level0) |>
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

phmap_cols_center <- RColorBrewer::brewer.pal(length(unique(metadata$center)), "Set2")
names(phmap_cols_center) <- unique(metadata[["center"]])

phmap_cols_axon_normal <- viridis::viridis(length(unique(metadata$log_axon_normal)))
names(phmap_cols_axon_normal) <- unique(metadata[["log_axon_normal"]])

phmap_list <-
    list(
        level0 = phmap_cols_level0,
        level2 = phmap_cols_level2,
        INCAT = phmap_cols_INCAT,
        center = phmap_cols_center,
        log_axon_normal = phmap_cols_axon_normal
    )

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

pdf(file.path("results", "mrvi", "heatmap_mrvi_average_all.pdf"), width = 10, height = 10)
print(phmap_mrvi)
dev.off()


# for each cluster separately ---

mrviPlot <- function(cluster) {
    mrvi_input <-
        read.csv(file.path("results", "mrvi", paste0("mrvi_cluster_", cluster, "_average.csv"))) |>
        tibble::column_to_rownames("X")


    metadata <-
        read_csv(file.path("lookup", "sample_lookup.csv")) |>
        select(level2, INCAT, center) |>
        mutate(level0 = if_else(level2 == "CTRL", "CTRL", "PNP")) |>
        data.frame()



    # make sure that rownames and columnnames of data and metadata are the same
    rownames(metadata) <- colnames(mrvi_input)

    # create colors for annotations
    phmap_cols_level2 <- pals::cols25(length(unique(metadata$level2)))
    names(phmap_cols_level2) <- unique(metadata[["level2"]])

    phmap_cols_level0 <- RColorBrewer::brewer.pal(length(unique(metadata$level0)), "Set1")[1:2]
    names(phmap_cols_level0) <- unique(metadata[["level0"]])

    phmap_cols_INCAT <- c("white", RColorBrewer::brewer.pal(length(unique(metadata$INCAT)) - 1, "Reds"))
    names(phmap_cols_INCAT) <- c("-", "1", "2", "3", "4", "5", "6")

    phmap_cols_center <- RColorBrewer::brewer.pal(length(unique(metadata$center)), "Set2")
    names(phmap_cols_center) <- unique(metadata[["center"]])

    phmap_list <-
        list(
            level0 = phmap_cols_level0,
            level2 = phmap_cols_level2,
            INCAT = phmap_cols_INCAT,
            center = phmap_cols_center
        )

    phmap_mrvi <-
        pheatmap(
            mrvi_input,
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
            annotation_row = metadata,
            annotation_colors = phmap_list,
            main = cluster,
        )


    pdf(file.path("results", "mrvi", paste0("heatmap_mrvi_cluster_", cluster, "_average.pdf")), width = 10, height = 10)
    print(phmap_mrvi)
    dev.off()
}

lapply(sc_merge@misc$cluster_order, mrviPlot)
