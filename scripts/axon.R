# analysis of axon counts

# libraries  ----
library(tidyverse)
library(patchwork)
library(Seurat)
library(BPCells)

axon_count_table <-
    readxl::read_xlsx(file.path("lookup", "axon_count_v2.xlsx")) |>
    dplyr::filter(is.na(remove)) |>
    dplyr::filter(sample != "NA") |>
    mutate(across(c(fascicle, normal_myelin:total_myelinated_axons, area_micrometer), parse_number)) |>
    mutate(fascicle = str_extract(fascicle, "[0-9]+")) 

axon_count_fascicle <- 
    axon_count_table |>
    group_by(sample, fascicle) |>
    summarize(across(c(normal_myelin:total_myelinated_axons, area_micrometer), sum), .groups = "drop") |>
    mutate(across(normal_myelin:total_myelinated_axons, function(x) x / area_micrometer, .names = "{.col}_density")) |>
    left_join(sample_lookup, by = "sample") |>
    mutate(level2 = factor(level2)) 

axon_count_sum <-
    axon_count_table |>
    group_by(sample) |>
    summarize(across(c(normal_myelin:total_myelinated_axons, area_micrometer), sum), .groups = "drop") |>
    mutate(across(normal_myelin:total_myelinated_axons, function(x) x / area_micrometer, .names = "{.col}_density")) |>
    left_join(sample_lookup, by = "sample") |>
    mutate(level2 = factor(level2))

axonPlot <- function(df, y_value) {
    df |>
        mutate(order = as.numeric(level2)) |>
        mutate(sample = fct_reorder(sample, order)) |>
        ggplot(aes(x = sample, y = .data[[y_value]], fill = level2, color = level2)) +
        geom_boxplot() +
        geom_jitter(height = 0, width = 0.1 ) +
        # scale_fill_manual(values = sc_merge@misc$level2_cols) +
        scale_color_manual(values = sc_merge@misc$level2_cols) +
        guides(fill = guide_legend(title = NULL)) +
        theme_classic() +
        ylab("") +
        xlab("") +
        ggtitle(y_value) +
        theme(axis.text.x = element_text(
            angle = 90,
            hjust = 1, vjust = 0.3
        ))
}

axon_plots_fascicle <- lapply(
    c(
        "normal_myelin",
        "thin_myelin",
        "very_small_axons",
        "total_myelinated_axons",
        "normal_myelin_density",
        "thin_myelin_density",
        "very_small_axons_density",
        "total_myelinated_axons_density"
    ),
    function(x) {
        axonPlot(df = axon_count_fascicle, y_value = x)
    }
)

axon_plots_fascicle_patch <- patchwork::wrap_plots(axon_plots_fascicle, ncol = 2)
ggsave(plot = axon_plots_fascicle_patch, file.path("results", "histo", "axon_count_fascicle.pdf"), width = 15, height = 15)

axon_plots_sum <- lapply(
    c(
        "normal_myelin",
        "thin_myelin",
        "very_small_axons",
        "total_myelinated_axons",
        "normal_myelin_density",
        "thin_myelin_density",
        "very_small_axons_density",
        "total_myelinated_axons_density"
    ),
    function(x) {
        axonPlot(df = axon_count_sum, y_value = x)
    }
)

axon_plots_sum_patch <- patchwork::wrap_plots(axon_plots_sum, ncol = 2)
ggsave(plot = axon_plots_sum_patch, file.path("results", "histo", "axon_count_sum.pdf"), width = 15, height = 15)

axonEphysioPlot <- function(x_axis, y_axis) {
    axon_count_fascicle |>
        group_by(sample) |>
        summarize(
            x_axis = median(.data[[x_axis]]),
            y_axis = median(.data[[y_axis]])
        ) |>
        ggplot(aes(x = x_axis, y = y_axis)) +
        geom_smooth(method = "lm") +
        geom_jitter(height = 0, width = 0.1) +
        guides(fill = guide_legend(title = NULL)) +
        theme_bw() +
        xlab(x_axis) +
        ylab(y_axis) 
}

x_axis_var <- axon_count_fascicle |>
    select(cmap_ulnar:ncv_sural) |>
    names()

# normal axons
normal_axon_ephysio_plots <- lapply(
    x_axis_var,
    function(x) {
        axonEphysioPlot(y_axis = "normal_myelin", x_axis = x)
    }
)

normal_axon_ephysio_plots_patch <- patchwork::wrap_plots(normal_axon_ephysio_plots, ncol = 3)
ggsave(plot = normal_axon_ephysio_plots_patch, file.path("results", "histo", "normal_axon_ephysio.pdf"), width = 15, height = 15)

# thinly myelinated axons
thin_axon_ephysio_plots <- lapply(
    x_axis_var,
    function(x) {
        axonEphysioPlot(y_axis = "thin_myelin", x_axis = x)
    }
)

thin_axon_ephysio_plots_patch <- patchwork::wrap_plots(thin_axon_ephysio_plots, ncol = 3)
ggsave(plot = thin_axon_ephysio_plots_patch, file.path("results", "histo", "thin_axon_ephysio.pdf"), width = 15, height = 15)


# normal myelin density
normal_density_axon_ephysio_plots <- lapply(
    x_axis_var,
    function(x) {
        axonEphysioPlot(y_axis = "normal_myelin_density", x_axis = x)
    }
)

normal_density_axon_ephysio_plots_patch <- patchwork::wrap_plots(normal_density_axon_ephysio_plots, ncol = 3)
ggsave(plot = normal_density_axon_ephysio_plots_patch, file.path("results", "histo", "normal_density_axon_ephysio.pdf"), width = 15, height = 15)

axon_plots_fascicle_ncv <- lapply(
    c(
        "normal_myelin",
        "thin_myelin",
        "very_small_axons",
        "total_myelinated_axons",
        "normal_myelin_density",
        "thin_myelin_density",
        "very_small_axons_density",
        "total_myelinated_axons_density"
    ),
    function(x) {
        axonPlotNCV(df = axon_count_fascicle, var = x)
    }
)

axon_plots_fascicle_ncv_patch <- patchwork::wrap_plots(axon_plots_fascicle_ncv, ncol = 2)
ggsave(plot = axon_plots_fascicle_ncv_patch, file.path("results", "histo", "axon_ncv_fascicle.pdf"), width = 15, height = 15)


axon_plots_sum_ncv <- lapply(
    c(
        "normal_myelin",
        "thin_myelin",
        "very_small_axons",
        "total_myelinated_axons",
        "normal_myelin_density",
        "thin_myelin_density",
        "very_small_axons_density",
        "total_myelinated_axons_density"
    ),
    function(x) {
        axonPlotNCV(df = axon_count_sum, var = x)
    }
)

axon_plots_sum_ncv_patch <- patchwork::wrap_plots(axon_plots_sum_ncv, ncol = 2)
ggsave(plot = axon_plots_sum_ncv_patch, file.path("results", "histo", "axon_ncv_sum.pdf"), width = 15, height = 15)

axon_count_mean <-
    axon_count_fascicle |>
    group_by(sample) |>
    dplyr::summarize(
        axon_normal = mean(normal_myelin),
    ) |>
    mutate(log_axon_normal = log(axon_normal))

sc_merge@meta.data <-
    sc_merge@meta.data |>
    tibble::rownames_to_column("barcode") |>
    dplyr::left_join(axon_count_mean, by = "sample") |>
    tibble::column_to_rownames(var = "barcode")

ic@meta.data <-
    ic@meta.data |>
    tibble::rownames_to_column("barcode") |>
    dplyr::left_join(axon_count_mean, by = "sample") |>
    tibble::column_to_rownames(var = "barcode")
