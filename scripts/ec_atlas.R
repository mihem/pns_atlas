library(Seurat)
library(tidyverse)
library(data.table)
library(qs)

# downloaded from here https://endotheliomics.shinyapps.io/ec_atlas/
counts_ec <- as.matrix(fread(file.path("raw_ec", "Data.csv")), rownames = 1)
counts_ec_matrix <- as(counts_ec, "dgCMatrix")
metadata_ec <- as.matrix(fread(file.path("raw_ec", "Metadata.csv")), rownames = 1)

# create seurat object and add metadata ---
ec_atlas <- CreateSeuratObject(counts = counts_ec_matrix, min.cells = 3, min.features = 200)
ec_atlas <- AddMetaData(ec_atlas, metadata_ec)

# run seurat pipeline ---
ec_atlas <- 
    ec_atlas  |>
    NormalizeData() |>
    ScaleData() |>
    FindVariableFeatures() |>
    RunPCA()

# better names
ec_atlas$cluster <- ec_atlas$Cluster
ec_atlas$endothelial <- ec_atlas$`Endothelial cell`
ec_atlas$tissue <- ec_atlas$Tissue

ec_atlas$Cluster <- NULL
ec_atlas$`Endothelial cell` <- NULL
ec_atlas$Tissue <- NULL

qs::qsave(ec_atlas, file.path("objects", "ec_atlas.qs"))




qs::qsave(ec_atlas, file.path("objects", "ec_atlas.qs"))

