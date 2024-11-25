#===============================================================================
# Xenium Data Analysis Pipeline
#
# Purpose: Process and analyze Xenium spatial transcriptomics data by:
# 1. Loading and preprocessing Xenium data
# 2. Integrating with single-nucleus data for cell type prediction
# 3. Generating spatial visualization plots
# 4. Analyzing marker abundance
#===============================================================================

# analyze xenium data

# load libraries
library(Seurat)
library(tidyverse)
library(scMisc)
library(sf)
future::plan("multicore", workers = 6)

# load Xenium data ----
xenium_paths <- list.dirs(file.path("xenium", "raw"), full.names = TRUE, recursive = FALSE)

xenium_meta <- 
    readxl::read_excel(file.path("lookup", "xenium_meta.xlsx"))  |>
    mutate(File = gsub(x = File, pattern = "(.+).tar", replacement = "\\1")) 

sample_lookup <- 
  readr::read_csv(file.path("lookup", "sample_lookup.csv")) |>
  janitor::clean_names() |>
  mutate(age_calc = lubridate::time_length(difftime(nerve_date, birth_date), "years")) |>
  mutate(age_calc = floor(age_calc)) |>
  mutate(age = coalesce(age_calc, age)) |>
  dplyr::select(-age_calc) |>
  mutate(level0 = if_else(level1 == "CTRL", "CTRL", "PNP"))

xenium_names <-
    tibble(File = basename(xenium_paths)) |>
    left_join(xenium_meta, join_by(File)) |>
    mutate(sample = str_extract(Name, "S\\d{2}")) |>
    select(sample) |>
    left_join(sample_lookup, join_by(sample))  |>
    mutate(name = paste0(sample, "_", level2)) |>
    pull(name)

xenium_objects <- lapply(
    xenium_paths,
    function(x) {
        LoadXenium(x, fov = "fov")
    }
) |>
    setNames(xenium_names)

# Load and preprocess Xenium data
# - Reads Xenium data from raw files
# - Filters cells with < 10 counts
# - Performs SCTransform normalization

# filter cells with less than 10 counts
xenium_objects <- lapply(xenium_objects, function(x) {
    x <- subset(x, subset = nCount_Xenium > 9)
})

#integration with seurat
for (xenium_object in (xenium_objects)) {
    xenium_object <- SCTransform(xenium_object, assay = "Xenium")
    xenium_object <- RunPCA(xenium_object)
}

sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)
sc_merge_small <- subset(sc_merge, downsample = 1000)

sc_merge_small <- SCTransform(sc_merge_small, verbose = FALSE)

predictions <- list()

# Integration with single-nucleus reference data
# - Finds transfer anchors between Xenium and snRNA-seq data
# - Transfers cell type labels from snRNA-seq to Xenium cells

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

# Visualization of spatial data
# - Generates spatial plots colored by predicted cell types
# - Creates cropped views of specific regions
# - Plots different cell populations (SC, T_NK cells)

for (i in names(xenium_objects)) {
    p1 <- ImageDimPlot(xenium_objects[[i]], group.by = "sn_predictions", cols = xenium_objects[[i]]@misc$cluster_col)
    ggsave(plot = p1, filename = file.path("results", "xenium", paste0("seurat_predict_", i, ".png")), width = 15, height = 15)
}

for (i in names(xenium_objects)) {
    Idents(xenium_objects[[i]]) <- xenium_objects[[i]]$sn_predictions
}

# predictions, sc only
for (i in c("S14_CIAP", "S30_VN")) {
    p1 <- ImageDimPlot(
        xenium_objects[[i]],
        cells = WhichCells(xenium_objects[[i]], idents = c("repairSC", "nmSC", "mySC")),
        group.by = "sn_predictions",
        cols = xenium_objects[[i]]@misc$cluster_col,
        axes = TRUE,
        dark.background = FALSE,
        size = 1
    ) + 
        theme_classic() + 
        theme(legend.title = element_blank())
    ggsave(plot = p1, filename = file.path("results", "xenium", "SC", paste0("seurat_predict_", i, ".png")), width = 8, height = 8, dpi = 1200)
}

for (i in c("S01_CIDP", "S24_CTRL")) {
    p1 <- ImageDimPlot(
        xenium_objects[[i]],
        cells = WhichCells(xenium_objects[[i]], idents = c("nmSC", "mySC")),
        group.by = "sn_predictions",
        cols = xenium_objects[[i]]@misc$cluster_col,
        axes = TRUE,
        dark.background = FALSE,
        size = 1
    ) + 
        theme_classic() + 
        theme(legend.title = element_blank())
    ggsave(plot = p1, filename = file.path("results", "xenium", "SC", paste0("seurat_predict_", i, ".png")), width = 8, height = 8, dpi = 1200)
}


# crop to representative areas
cropped_coords_S24_CTRL <- Crop(xenium_objects[["S24_CTRL"]][["fov"]], x = c(1900, 2800), y = c(2000, 3200), coords = "plot")
xenium_objects[["S24_CTRL"]][["zoom"]] <- cropped_coords_S24_CTRL

cropped_coords_S30_VN <- Crop(xenium_objects[["S30_VN"]][["fov"]], x = c(500, 1800), y = c(2000, 3200), coords = "plot")
xenium_objects[["S30_VN"]][["zoom"]] <- cropped_coords_S30_VN

for (i in c("S30_VN")) {
    p1 <- ImageDimPlot(
        xenium_objects[[i]],
        cells = WhichCells(xenium_objects[[i]], idents = c("repairSC", "nmSC", "mySC")),
        group.by = "sn_predictions",
        cols = xenium_objects[[i]]@misc$cluster_col,
        # axes = TRUE,
        fov = "zoom",
        size = 1, 
        flip_xy = FALSE) + 
        theme(legend.title = element_blank())
    ggsave(plot = p1, filename = file.path("results", "xenium", "SC", paste0("seurat_predict_zoom_", i, ".png")), width = 3.5, height = 3.5, dpi = 600)
}

for (i in c("S24_CTRL")) {
    p1 <- ImageDimPlot(
        xenium_objects[[i]],
        cells = WhichCells(xenium_objects[[i]], idents = c("nmSC", "mySC")),
        group.by = "sn_predictions",
        cols = xenium_objects[[i]]@misc$cluster_col,
        # axes = TRUE,
        fov = "zoom",
        size = 1, 
        flip_xy = FALSE) + 
        theme(legend.title = element_blank())
    ggsave(plot = p1, filename = file.path("results", "xenium", "SC", paste0("seurat_predict_zoom_", i, ".png")), width = 3.5, height = 3.5, dpi = 600)
}


# only T_NK
for (i in names(xenium_objects)) {
    try({
        p1 <- ImageDimPlot(
            xenium_objects[[i]],
            cells = WhichCells(xenium_objects[[i]], idents = c("T_NK")),
            group.by = "sn_predictions",
            cols = pals::cols25()[2],
            size = 1.0,
            dark.background = FALSE
        ) +
            theme_classic() +
            theme(legend.title = element_blank());
        ggsave(plot = p1, filename = file.path("results", "xenium", "t_nk", paste0("seurat_predict_", i, ".png")), width = 8, height = 8, dpi = 1200)
    })
}

# Analysis of cell type abundances and marker quantification
# - Summarizes predictions into broader cell type groups
# - Analyzes manually quantified marker abundance
# - Generates statistical plots and comparisons

# summarize predictions to larger groups
predictions_lookup <-
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

# S24_CTRL cropped 
cropped_coords_S24_CTRL_group <- Crop(xenium_objects[["S24_CTRL"]][["fov"]], x = c(1500, 2500), y = c(2200, 3500), coords = "plot")
xenium_objects[["S24_CTRL"]][["zoom"]] <- cropped_coords_S24_CTRL_group
S24_group_plot <- ImageDimPlot(xenium_objects[["S24_CTRL"]], group.by = "sn_predictions_group", cols = prediction_group_col, axes = FALSE, size = 1, fov = "zoom") + 
    theme(legend.title = element_blank())
ggsave(plot = S24_group_plot, filename = file.path("results", "xenium", paste0("seurat_predict_group_S24_CTRL_cropped.png")), width = 3.5, height = 3.5, dpi = 600)


for (i in names(xenium_objects)) {
    Idents(xenium_objects[[i]]) <- xenium_objects[[i]]$sn_predictions_group
}

qs::qsave(xenium_objects, file.path("objects", "xenium_objects.qs"))
xenium_objects <- qs::qread(file.path("objects", "xenium_objects.qs"))


# abundance of xenium markers based on manual quantification ----
quanti <-
    read_csv(file.path("lookup", "xenium_manual_quantification.csv"))|>
    pivot_longer(!sample, names_to = "name", values_to = "value") |>
    separate_wider_delim(name, delim = "_", names = c("variable", "area"))

quanti_stats <-
    quanti |>
    dplyr::filter(area != "all") |>
    group_by(sample, variable) |>
    summarize(
        mean = mean(value, na.rm = TRUE),
        sum = sum(value, na.rm = TRUE)
    ) |>
    pivot_longer(c(mean, sum), names_to = "name", values_to = "value") |>
    unite("variable", variable, name)


quanti_result <-
    quanti |>
    dplyr::filter(area == "all") |>
    select(sample, variable, value) |>
    bind_rows(quanti_stats)  |>
    arrange(sample) |>
    pivot_wider(names_from = variable, values_from = value) |>
    left_join(select(sample_lookup, sample, level2), join_by(sample)) |>
    mutate(level2 = factor(level2, levels = sc_merge@misc$level2_order)) |>
    mutate(epiarea = nervearea - endoarea_sum) |>
    mutate(epiPTPRC = nervePTPRC - endoPTPRC_sum) |>
    mutate(epiMS4A1 = nerveMS4A1 - endoMS4A1_sum) |>
    mutate(epiCD3E = nerveCD3E - endoCD3E_sum) |>
    mutate(epiPTPRC_density = epiPTPRC / epiarea) |>
    mutate(epiMS4A1_density = epiMS4A1 / epiarea) |>
    mutate(epiCD3E_density = epiCD3E / epiarea) |>
    mutate(endoPTPRC_density_sum = endoPTPRC_sum / endoarea_sum) |>
    mutate(endoMS4A1_density_sum = endoMS4A1_sum / endoarea_sum) |>
    mutate(endoCD3E_density_sum = endoCD3E_sum / endoarea_sum)

plotQuanti <- function(y_value, name) {
    quanti_result |>
        ggplot(aes(x = level2, y = .data[[y_value]], fill = level2)) +
        geom_boxplot() +
        geom_jitter(height = 0, width = 0.1) +
        scale_fill_manual(values = sc_merge@misc$level2_cols) +
        guides(fill = guide_legend(title = NULL)) +
        theme_classic() +
        ylab("") +
        xlab("") +
        ggtitle(name) +
        theme(
            legend.position = "none",
            axis.text.x = element_text(
                angle = 90,
                hjust = 1,
                vjust = 0.3,
            )
        )
    ggsave(file.path("results", "xenium", "quantification", paste0(name, ".pdf")), width = 1.75, height = 2.5)
}

plotQuanti(y_value = "epiPTPRC_density", name = "PTPRC_epi_density")
plotQuanti(y_value = "endoPTPRC_density_sum", name = "PTPRC_endo_density_sum")


plotQuanti(y_value = "epiMS4A1_density", name = "MS4A1_epi_density")
plotQuanti(y_value = "endoMS4A1_density_sum", name = "MS4A1_endo_density_sum")

plotQuanti(y_value = "epiCD3E_density", name = "CD3E_epi_density")
plotQuanti(y_value = "endoCD3E_density_sum", name = "CD3E_endo_density_sum")

# abundance of xenium cell predictions ----
cells_list <-
    lapply(xenium_objects, function(x) (x$sn_predictions)) |>
    unlist()

cells_predicted <-
    tibble(cluster =  unname(cells_list), sample = names(cells_list)) |>
    mutate(condition = str_replace(sample, "(S\\d+)_(\\w+).*", "\\2")) |>
    mutate(sample = str_replace(sample, "(S\\d+)_.*", "\\1"))  |>
    dplyr::filter(!is.na(cluster))

# boxplot sc
boxplot_xenium_sc_t_nk <-
    cells_predicted |>
    dplyr::count(cluster, sample) |>
    pivot_wider(names_from = sample, values_from = n) |>
    mutate(across(where(is.numeric), function(x) x / sum(x, na.rm = TRUE) * 100)) |>
    pivot_longer(!cluster, names_to = "sample", values_to = "percent") |>
    left_join(select(sample_lookup, sample, level2)) |>
    rename(condition = level2) |>
    dplyr::filter(cluster %in% c("mySC", "nmSC", "repairSC", "T_NK"))  |>
    mutate(percent  = replace_na(percent, 0)) |>
    mutate(condition = factor(condition, levels = sc_merge@misc$level2_order)) |>
    ggplot(aes(x = condition, y = percent, fill = condition)) +
    geom_boxplot() +
    geom_point() + 
    theme_classic() +
    facet_wrap(vars(cluster), scales = "free_y", nrow = 1) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
    xlab("") +
    ylab("percentage") +
    scale_fill_manual(values = sc_merge@misc$level2_cols) +
    theme(legend.position = "none")


ggsave(file.path("results", "xenium", "abundance_xenium_sc_t_nk.pdf"),
    plot = boxplot_xenium_sc_t_nk,
    width = 5,
    height = 3
)


# Export results ----
# export prediction annotation for python analysis
for (i in names(xenium_objects)) {
    arrow::write_parquet(xenium_objects[[i]]@meta.data, sink = file.path("results", "xenium", paste0("xenium_predictions_", i, ".parquet")))
}
