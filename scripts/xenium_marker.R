# libraries  ----
library(Seurat)
library(BPCells)
library(SeuratObject)
library(tidyverse)
library(writexl)
library(patchwork)
library(conflicted)
library(qs)
library(pals)
library(scMisc)
library(Polychrome)
library(readxl)
library(DropletUtils)

# devtools::install_github("mihem/scMisc")
# detach(package:scMisc, unload = TRUE)

# general settings  ----
conflicts_prefer(base::setdiff)
my_cols_25 <- pals::cols25()
my_cols_50 <- unname(Polychrome::createPalette(50, pals::cols25()))

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)
sc_merge_small <- qs::qread(file.path("objects", "sc_merge_small.qs"))
sc_xenium <- subset(sc_merge, sample %in% c("S22", "S24", "S29", "S30"))

dplyr::count(sc_xenium@meta.data, cluster) |>
  dplyr::arrange(n)


markers_xenium <- read_csv(file.path("lookup", "xenium_list_jolien.csv"))

# create h5 file as a reference dataset
sc_merge_small$RNA$counts <- as(object = sc_merge_small[["RNA"]]$counts, Class = "dgCMatrix")
sc_xenium$RNA$counts <- as(object = sc_xenium[["RNA"]]$counts, Class = "dgCMatrix")

DropletUtils::write10xCounts(
  x = sc_merge_small$RNA$counts,
  path = file.path("objects", "sc_merge_small_count.h5"),
  version = "3"
)

DropletUtils::write10xCounts(
  x = sc_xenium$RNA$counts,
  path = file.path("objects", "sc_xenium_count.h5"),
  version = "3"
)

# DropletUtils::write10xCounts(
#   x = sc_merge_small$RNA$counts,
#   path = file.path("objects", "sc_merge_small_count"),
#   version = "3"
# )

# save annotations for Xenium
str(sc_merge_small@meta.data)

data.frame(annotation = sc_merge_small$cluster) |>
  rownames_to_column(var = "barcode") |>
  write_csv(file.path("objects", "sc_merge_small_annotation.csv"))

data.frame(annotation = sc_xenium$cluster) |>
  rownames_to_column(var = "barcode") |>
  write_csv(file.path("objects", "sc_xenium_annotation.csv"))

zip(
  zipfile = file.path("objects", "sc_xenium.zip"),
  files = c(file.path("objects", "sc_xenium_count.h5"),
            file.path("objects", "sc_xenium_annotation.csv"))
)


# dotplot
DefaultAssay(sc_merge) <- "RNA"

rownames(sc_merge_small)

lapply(names(markers_xenium)[-11],
       function(x) {
         dotPlot(
           path = file.path("lookup", "xenium_list_jolien.csv"),
           object = sc_merge,
           par = x,
           dot_min = 0.01,
           height = 10,
           width = 10
         )
       })
    
dotPlot(
  path = file.path("lookup", "xenium_list_jolien.csv"),
  object = sc_merge,
  par = "top_10_EC_Macro_SC_tactocyte",
  dot_min = 0.01,
  height = 10,
  width = 30
)

dotPlot(
  path = file.path("lookup", "markers.csv"),
  object = sc_xenium,
  par = "xenium_cluster",
  dot_min = 0.01,
  height = 7,
  width = 20
)

str(sc_xenium@meta.data)

Idents(sc_xenium) <- sc_xenium$level2

cluster_order_de <-
  expand_grid(
    cluster = sc_xenium@misc$cluster_order,
    level2 = unique(sc_xenium$level2)
  ) |>
  mutate(cluster_level2 = paste0(cluster, "_", level2))


Idents(sc_xenium) <- factor(paste0(sc_xenium$cluster, "_", sc_xenium$level2),
  levels = cluster_order_de$cluster_level2
)
                           

dotPlot(
  path = file.path("lookup", "markers.csv"),
  object = sc_xenium,
  par = "xenium_de",
  dot_min = 0.01,
  height = 10,
  width = 12
)

