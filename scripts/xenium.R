# analyse xenium data

# load libraries
library(Seurat)
library(tidyverse)

future::plan("multicore", workers = 6)

# load Xenium data

xenium_path <- "xenium/raw/output-XETG00051__LAB5319-Slide1__S154297__20240202__141517"
xenium_paths <- list.dirs(file.path("xenium", "raw"), full.names = TRUE, recursive = FALSE)

xenium_meta <- 
    readxl::read_excel(file.path("lookup", "xenium_meta.xlsx"))  |>
    mutate(File = gsub(x = File, pattern = "(.+).tar", replacement = "\\1")) 

xenium_names <-
    tibble(File = basename(xenium_paths)) |>
    left_join(xenium_meta) |>
    mutate(Name = str_extract(Name, "S\\d{2}")) |>
    pull(Name)

xenium_objects <- lapply(
    xenium_paths,
    function(x) {
        LoadXenium(x, fov = "fov")
    }
) |>
    setNames(xenium_names)

# function to plot Xenium image plots for all samples
XeniumImagePlot <- function(xenium_objects, markers_xenium, sel_cluster) {
    molecules_sel <-
        markers_xenium |>
        dplyr::filter(cluster %in% {{ sel_cluster }}) |>
        pull(transcript)
    for (i in seq_along(xenium_objects)) {
        dir.create(file.path("results", "xenium", names(xenium_objects)[i]))
        plot <- ImageDimPlot(xenium_objects[[i]],
            fov = "fov",
            molecules = molecules_sel,
            nmols = 20000,
            cols = "white"
        )
        ggsave(
            plot = plot,
            filename = file.path("results", "xenium", names(xenium_objects)[i], paste0(sel_cluster, ".png")),
            width = 8,
            height = 8
        )
    }
}

markers_xenium <- read_xlsx(file.path("lookup", "xenium_list_final.xlsx"))

clusters_xenium <- unique(markers_xenium$cluster)

lapply(
    clusters_xenium,
    function(x) {
        XeniumImagePlot(xenium_objects = xenium_objects, markers_xenium = markers_xenium, sel_cluster = x)
    })


# deconvolution using spacexr
remotes::install_github("dmcable/spacexr")

library(spacexr)

xenium.obj <- xenium_objects[[4]]


query.counts <- GetAssayData(xenium.obj, assay = "Xenium", slot = "counts")
coords <- GetTissueCoordinates(xenium.obj, which = "centroids")
rownames(coords) <- coords$cell
coords$cell <- NULL
query <- SpatialRNA(coords, query.counts, colSums(query.counts))


Idents(sc_xenium) <- sc_xenium$cluster
sc_xenium_small <- subset(sc_xenium, downsample = 1000)
#remove macro1 becuase less than 25 cells
sc_xenium_small <- subset(sc_xenium_small, cluster %in% c("Macro1"), invert = TRUE)
sc_xenium_small$cluster <- droplevels(sc_xenium_small$cluster)

counts <- GetAssayData(sc_xenium_small, assay = "RNA", slot = "counts")

cluster <- sc_xenium_small$cluster
nUMI <- sc_xenium_small$nCount_RNA

reference <- Reference(counts, cluster, nUMI)

# run RCTD with many cores (~10 min)
RCTD <- create.RCTD(query, reference, max_cores = 4)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

scMisc::lss()

annotations.df <- RCTD@results$results_df

dplyr::count(annotations.df, first_type, sort = TRUE)
annotations <- annotations.df$first_type
names(annotations) <- rownames(annotations.df)
xenium.obj$predicted.celltype <- annotations

dplyr::count(xenium.obj@meta.data, predicted.celltype, sort = TRUE)
qs::qsave(xenium.obj, file.path("objects", "xenium.qs")

p1 <- ImageDimPlot(xenium.obj, group.by = "predicted.celltype", cols = "polychrome")
ggsave(plot = p1, filename = "test.png")

# niche assay
xenium.obj <- BuildNicheAssay(
    object = xenium.obj,
    fov = "fov",
    # group.by = "predicted.celltype",
    niches.k = 5,
    neighbors.k = 30
)

celltype.plot <- ImageDimPlot(xenium.obj, group.by = "predicted.celltype", size = 1.5, cols = "polychrome",
    dark.background = F) + ggtitle("Cell type")
niche.plot <- ImageDimPlot(xenium.obj, group.by = "niches", size = 1.5, dark.background = F) + ggtitle("Niches") +
    scale_fill_manual(values = c("#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B"))
celltype.plot | niche.plot