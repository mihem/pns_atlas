# analyse xenium with SPIA

# load libraries
library(Seurat)
library(tidyverse)
library(scMisc)
library(SPIAT)
library(patchwork)
library(RColorBrewer)
library(pheatmap)
future::plan("multicore", workers = 6)

# load preprocessed data
xenium_objects <- qs::qread(file.path("objects", "xenium_objects.qs"))

# create SpatialExperiment object from Seurat ----
seurat_to_spe <- function(seurat_object) {
    seurat_object$sn_predictions <- as.character(seurat_object$sn_predictions)
    seurat_object$sn_predictions[is.na(seurat_object$sn_predictions)] <- "unknown"
    seurat_object$sn_predictions_group <- as.character(seurat_object$sn_predictions_group)
    seurat_object$sn_predictions_group[is.na(seurat_object$sn_predictions_group)] <- "unknown"
    intensity_matrix <- as.matrix(GetAssayData(seurat_object, assay = "Xenium", layer = "counts"))
    coords <- GetTissueCoordinates(seurat_object, which = "centroids")
    phenotypes <- seurat_object$sn_predictions
    spe_object <- format_image_to_spe(
        format = "general",
        intensity_matrix = intensity_matrix,
        phenotype = phenotypes,
        coord_x = coords$x,
        coord_y = coords$y
    )
    spe_object$prediction_group <- seurat_object$sn_predictions_group
    return(spe_object)
}

spe_objects <- lapply(xenium_objects, seurat_to_spe)

spiat_cell_categories <-
    lapply(names(spe_objects),
        FUN = function(name) {
            plot_cell_categories(
                spe_object = spe_border[[name]],
                # categories_of_interest = c("endoC"),
                categories_of_interest = c("SC", "endoC", "periC", "epiC", "VSMC", "PC", "EC", "IC"),
                colour_vector = prediction_group_col,
                feature_colname = "prediction_group",
                cex = 0.1
            ) +
            ggtitle(name)
        }
    )

spiat_cell_catgories_plot <- wrap_plots(spiat_cell_categories, ncol = 1)

ggsave(file.path("results", "xenium", "spiat", "cell_categories.pdf"),
    plot = spiat_cell_catgories_plot,
    width = 7, height = 30
)

# # distance between cell types
# # takes some times to calculcate
# distances <-
#     calculate_pairwise_distances_between_celltypes(
#         spe_object = general_format_image,
#         cell_types_of_interest = c("epiC", "periC", "endoC", "SC"),
#         feature_colname = "prediction_group"
#     )

# plot_cell_distances_violin(distances)
# summary_distances <- calculate_summary_distances_between_celltypes(distances)

# p1 <- plot_distance_heatmap(phenotype_distances_result = summary_distances, metric = "mean")
# ggsave(file.path("results", "xenium", "spiat", "distance_heatmap_s01.pdf"), plot = p1, width = 4, height = 3)


# # mixing scores for cellular interactions
# mixing_score_sc <- list()
# for (i in c("endoC", "periC", "epiC")) {
#     mixing_score_sc[[i]] <-
#         mixing_score_summary(
#             spe_object = general_format_image,
#             reference_celltype = "SC",
#             target_celltype = i,
#             radius = 30,
#             feature_colname = "prediction_group"
#         )
# }

# p1 <- mixing_score_sc |>
#     bind_rows() |>
#     ggplot(aes(x = Target, y = Normalised_mixing_score, fill = Target)) +
#     geom_col() + 
#     theme_classic() +
#     scale_fill_manual(values = prediction_group_col) +
#     ggtitle("Mixing score to SC") + 
#     theme(legend.position = "none")

# ggsave(file.path("results", "xenium", "spiat", "mixing_score_sc_s01.pdf"),
#     plot = p1, width = 3, height = 3
# )

# spatial heterogeneity
# test1 <-
#     calculate_entropy(
#         general_format_image,
#         cell_types_of_interest = c("SC", "epiC", "periC", "endoC"),
#         feature_colname = "prediction_group"
#     )

gridEntropy <- function(name) {
    grid <- grid_metrics(
        spe_objects[[name]],
        FUN = calculate_entropy,
        n_split = 30,
        cell_types_of_interest = c("periC", "endoC"),
        feature_colname = "prediction_group"
    )
    pdf(file.path("results", "xenium", "spiat", paste0(name, "_grid_entropy.pdf")),
        width = 8, height = 8
    )
    plot(grid)
    dev.off()
    return(grid)
}

spe_grid <- lapply(names(spe_objects), gridEntropy)

grid_entropy_pct <- lapply(spe_grid, function(x) calculate_percentage_of_grids(x, threshold = 0.75, above = TRUE))
grid_entropy_autocorrelation <- lapply(spe_grid, function(x) calculate_spatial_autocorrelation(x, metric = "globalmoran"))

grid_entropy_results <- 
    tibble(
        sample = names(spe_objects),
        grid_entropy_pct = map_dbl(spe_grid, function(x) calculate_percentage_of_grids(x, threshold = 0.75, above = TRUE)),
        grid_entropy_autocorrelation = map_dbl(spe_grid, function(x) calculate_spatial_autocorrelation(x, metric = "globalmoran"))
    )

grid_entropy_pct_plot <- 
    grid_entropy_results |>
    ggplot(aes(x = sample, y = grid_entropy_pct, fill = sample)) +
    geom_col() +
    theme_classic() + 
    theme(legend.position = "none") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    scale_fill_manual(values = pals::cols25())

ggsave(file.path("results", "xenium", "spiat", "grid_entropy_pct.pdf"),
    plot = grid_entropy_pct_plot, width = 2, height = 3
)

grid_entropy_autocorrelation <-
    grid_entropy_results |>
    ggplot(aes(x = sample, y = grid_entropy_autocorrelation, fill = sample)) +
    geom_col() +
    theme_classic() +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    scale_fill_manual(values = pals::cols25())

ggsave(file.path("results", "xenium", "spiat", "grid_entropy_autocorrelation.pdf"),
    plot = grid_entropy_autocorrelation, width = 2, height = 3
)

spe_border <- lapply(
    spe_objects,
    function(x) {
        identify_bordering_cells(
            x,
            reference_cell = "endoC",
            feature_colname = "prediction_group",
        )
    }
)

phmap_colors <- colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100)
names(spe_border)

spia_hmap_border <- function(name) {
    phmap <-
        table(spe_border[[name]]@colData$prediction_group, spe_border[[name]]@colData$Region) |>
        as.data.frame.matrix() |>
        data.frame() |>
        rownames_to_column("cluster") |>
        dplyr::filter(cluster != "unknown") |>
        column_to_rownames("cluster") |>
        pheatmap(
            cluster_cols = FALSE,
            scale = "row",
            color = phmap_colors,
            cellwidth = 10,
            cellheight = 10,
            treeheight_row = 15,
            treeheight_col = 15,
            clustering_distance_rows = "euclidean",
            clustering_method = "ward.D2",
            border_color = NA,
            cutree_rows = 2,
            main = name
        )
        pdf(file.path("results", "xenium", "spiat", paste0("phmap_border_", name, ".pdf")),
            width = 3, height = 5
        )
        print(phmap)
        dev.off()
}

lapply(names(spe_border), spia_hmap_border) 

# s01_distance <- calculate_distance_to_margin(s01_border)

# s01_structure <-
#     define_structure(
#         s01_distance,
#         cell_types_of_interest = c("epiC", "periC", "endoC", "IC", "SC"),
#         feature_colname = "prediction_group",
#         n_margin_layers = 5
#     )

# plot_cell_categories(s01_structure, feature_colname = "Structure")
# ggsave(file.path("results", "xenium", "spiat", "s01_structure.pdf"),)

# colData(s01_structure)

# phmap_structure_s01 <-
#     table(s01_structure@colData$prediction_group, s01_structure@colData$Structure) |>
#     as.data.frame.matrix() |>
#     data.frame() |>
#     rownames_to_column("cluster") |>
#     dplyr::filter(cluster != "unknown") |>
#     column_to_rownames("cluster") |>
#     pheatmap(
#         cluster_cols = FALSE,
#         scale = "row",
#         color = phmap_colors,
#         cellwidth = 10,
#         cellheight = 10,
#         treeheight_row = 15,
#         treeheight_col = 15,
#         clustering_distance_rows = "euclidean",
#         clustering_method = "ward.D2",
#         border_color = NA,
#         cutree_rows = 4,
#     )

# pdf(file.path("results", "xenium", "spiat", "phmap_structure_s01.pdf"),
#     width = 3, height = 5
# )
# print(phmap_structure_s01)
# dev.off()
