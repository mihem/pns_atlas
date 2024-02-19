# analyse xenium data

# load libraries
library(Seurat)
library(tidyverse)

future::plan("multicore", workers = 6)

# load Xenium data

xenium_paths <- list.dirs(file.path("xenium", "raw"), full.names = TRUE, recursive = FALSE)

# Get the center position of each centroid. There is one row per cell in this dataframe.
str(GetTissueCoordinates(xenium_s11[["fov"]][["centroids"]]))

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


#integration with seurat
for (xenium_object in (xenium_objects)) {
    xenium_object <- subset(xenium_object, subset = nCount_Xenium > 0)
    xenium_object <- SCTransform(xenium_object, assay = "Xenium")
    xenium_object <- RunPCA(xenium_object)
}

sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)
sc_merge_small <- subset(sc_merge, downsample = 1000)

sc_merge_small <- SCTransform(sc_merge_small, verbose = FALSE)

predictions <- list()

for (i in names(xenium_objects)) {
    anchor <- FindTransferAnchors(
        reference = sc_merge_small,
        query = xenium_objects[[i]],
        normalization.method = "SCT"
    )
    predictions[[i]] <- TransferData(anchor, refdata = sc_merge_small$cluster)
}


#function to store predictions in seurat object
storePred <- function(predictions, label_col, score_col, seu_obj) {
  predictions_prep <-
    predictions |>
    tibble::rownames_to_column("barcode") |>
    dplyr::select(predicted.id, prediction.score.max, barcode) |>
    dplyr::mutate(predicted.id = ifelse(prediction.score.max < 0.3, "unknown", predicted.id)) |>
    tibble::as_tibble() |>
    dplyr::rename(
      {{ label_col }} := predicted.id,
      {{ score_col }} := prediction.score.max
    )

  seu_obj@meta.data <-
    seu_obj@meta.data |>
    tibble::rownames_to_column("barcode") |>
    dplyr::left_join(predictions_prep, by = "barcode") |>
    tibble::column_to_rownames(var = "barcode")

  return(seu_obj)
}

for (i in names(xenium_objects)) {
    xenium_objects[[i]] <- storePred(
        predictions = predictions[[i]],
        label_col = "sn_predictions",
        score_col = "sn_predictions_score",
        seu_obj = xenium_objects[[i]]
    )
}

for (i in names(xenium_objects)) {
    xenium_objects[[i]]@misc$cluster_order <- sc_merge@misc$cluster_order
    xenium_objects[[i]]@misc$cluster_col <- sc_merge@misc$cluster_col
    xenium_objects[[i]]$sn_predictions <- factor(xenium_objects[[i]]$sn_predictions, levels = xenium_objects[[i]]@misc$cluster_order)
}


for (i in names(xenium_objects)) {
    p1 <- ImageDimPlot(xenium_objects[[i]], group.by = "sn_predictions", cols = xenium_objects[[i]]@misc$cluster_col)
    ggsave(plot = p1, filename = file.path("results", "xenium", paste0("seurat_predict_", i, ".png")), width = 15, height = 15)
}


# summarize predictions to larger groups
predictions_lookup1 <-
    tibble(
        sn_predictions = xenium_objects[[1]]@misc$cluster_order,
        sn_predictions_group = c(rep("SC", 4), NA, "endoC", rep("periC", 3), "epiC", "VSMC", rep("PC", 2), rep("EC", 5), rep("IC", 6))
    )

for (i in names(xenium_objects)) {
    xenium_objects[[i]]@meta.data <-
        xenium_objects[[i]]@meta.data |>
        tibble::rownames_to_column("barcode") |>
        dplyr::left_join(predictions_lookup, by = "sn_predictions") |>
        tibble::column_to_rownames(var = "barcode")
}


prediction_group_level <- c("SC", "endoC", "periC", "epiC", "VSMC", "PC", "EC", "IC")
prediction_group_col <- setNames(pals::cols25(length(prediction_group_level)), prediction_group_level)

for (i in names(xenium_objects)) {
    xenium_objects[[i]]$sn_predictions_group <- factor(xenium_objects[[i]]$sn_predictions_group, levels = prediction_group_level)
}

for (i in names(xenium_objects)) {
    p1 <- ImageDimPlot(xenium_objects[[i]], group.by = "sn_predictions_group", cols = prediction_group_col)
    ggsave(plot = p1, filename = file.path("results", "xenium", paste0("seurat_predict_group_", i, ".png")), width = 15, height = 15)
}

qs::qsave(xenium_objects, file.path("objects", "xenium_objects.qs"))

# # deconvolution using spacexr
## comment: did not work as well as Integration with seurat, only few cells mapped
# xenium_s11 <- LoadXenium(xenium_paths[4])
# remotes::install_github("dmcable/spacexr")

# library(spacexr)



# query.counts <- GetAssayData(xenium_s11, assay = "Xenium", slot = "counts")
# coords <- GetTissueCoordinates(xenium_s11, which = "centroids")
# rownames(coords) <- coords$cell
# coords$cell <- NULL
# query <- SpatialRNA(coords, query.counts, colSums(query.counts))

# sc_merge_small <- subset(sc_merge, downsample = 1000)
# table(sc_merge_small$cluster)

# Idents(sc_xenium) <- sc_xenium$cluster
# sc_xenium_small <- subset(sc_xenium, downsample = 1000)
# #remove macro1 becuase less than 25 cells
# sc_xenium_small <- subset(sc_xenium_small, cluster %in% c("Macro1"), invert = TRUE)
# sc_xenium_small$cluster <- droplevels(sc_xenium_small$cluster)


# s11 <- subset(sc_xenium, subset = sample %in% c("S11"))
# #remove small clusters
# s11 <- subset(s11, cluster %in% c("damageSC", "periC3", "VSMC", "venEC", "Macro1", "Granulo", "B"), invert = TRUE)
# s11$cluster <- droplevels(s11$cluster)

# seurat_obj <- s11
# seurat_obj <- sc_merge_small

# dplyr::count(seurat_obj@meta.data, cluster, sort = TRUE)

# counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
# cluster <- seurat_obj$cluster
# nUMI <- seurat_obj$nCount_RNA
# reference <- Reference(counts, cluster, nUMI)

# # run RCTD with many cores (~10 min)
# RCTD <- create.RCTD(query, reference, max_cores = 4)
# RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

# annotations.df <- RCTD@results$results_df

# dplyr::count(annotations.df, first_type, sort = TRUE)
# annotations <- annotations.df$first_type
# names(annotations) <- rownames(annotations.df)
# xenium.obj$predicted.celltype <- annotations

# dplyr::count(xenium.obj@meta.data, predicted.celltype, sort = TRUE)
# qs::qsave(xenium.obj, file.path("objects", "xenium_s11_rctd.qs"))
# xenium.obj <- qs::qread(file.path("objects", "xenium_s11_rctd.qs"))

# keep.cells <- Cells(xenium.obj)[!is.na(xenium.obj$predicted.celltype)]
# xenium.obj <- subset(xenium.obj, cells = keep.cells)

# p1 <- ImageDimPlot(xenium.obj, group.by = "predicted.celltype", cols = "polychrome")
# ggsave(plot = p1, filename = "test.png")

# # does not work, probably too few cells predicted
# # # niche assay
# # xenium.obj <- BuildNicheAssay(
# #     object = xenium.obj,
# #     fov = "fov",
# #     # group.by = "predicted.celltype",
# #     niches.k = 5,
# #     neighbors.k = 30
# # )

# # celltype.plot <- ImageDimPlot(xenium.obj, group.by = "predicted.celltype", size = 1.5, cols = "polychrome",
# #     dark.background = F) + ggtitle("Cell type")
# # niche.plot <- ImageDimPlot(xenium.obj, group.by = "niches", size = 1.5, dark.background = F) + ggtitle("Niches") +
# #     scale_fill_manual(values = c("#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B"))
# # celltype.plot | niche.plot


# # find transfer anchors

