# analysis of axon counts

# libraries  ----
library(tidyverse)
library(patchwork)
library(Seurat)
library(BPCells)
library(qs)

# read data ----
sc_merge <- qread(file.path("objects", "sc_merge.qs"), nthread = 4)

sample_lookup <-
    read_csv(file.path("lookup", "sample_lookup.csv")) |>
    janitor::clean_names() |>
    dplyr::rename(
        ncv_tibial_motoric = ncv_tibial_motoric_in_m_s,
        cmap_tibial_motoric = cmap_tibial_in_m_v,
        f_latency_tibial = min_f_latency_tibial_in_ms,
        ncv_peroneal_motoric = ncv_peroneal_motoric_in_m_s,
        cmap_peroneal_motoric = cmap_peroneal_in_m_v,
        ncv_ulnar_motoric = ncv_ulnar_motoric_in_m_s,
        cmap_ulnar = cmap_ulnar_in_m_v,
        f_latency_ulnar = min_f_latency_ulnar_in_ms,
        snap_sural = snap_sural_in_m_v,
        ncv_sural = ncv_sural_in_m_s
    ) |>
    mutate(level0 = if_else(level1 == "CTRL", "CTRL", "PNP")) |>
    select(sample, level0, level2, incat, center, cmap_ulnar:ncv_sural) |>
    mutate(across(cmap_ulnar:ncv_sural, as.numeric))

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
    mutate(level2 = factor(level2, levels = sc_merge@misc$level2_order))

axon_count_sum <-
    axon_count_table |>
    group_by(sample) |>
    summarize(across(c(normal_myelin:total_myelinated_axons, area_micrometer), sum), .groups = "drop") |>
    mutate(across(normal_myelin:total_myelinated_axons, function(x) x / area_micrometer, .names = "{.col}_density")) |>
    left_join(sample_lookup, by = "sample") |>
    mutate(level2 = factor(level2))

# add normal axon counts to seurat objects ----
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

# plot normal axon counts grouped by level2 ----
axon_count_level2 <-
    axon_count_fascicle |>
    group_by(sample) |>
    mutate(axon_count = mean(normal_myelin)) |>
    select(sample, level2, axon_count) |>
    distinct() |>
    ggplot(aes(x = level2, y = axon_count, fill = level2)) +
    geom_boxplot() +
    geom_jitter(height = 0, width = 0.1) +
    scale_fill_manual(values = sc_merge@misc$level2_cols) +
    theme_classic() +
    ylab("") +
    xlab("") +
    ggtitle("normal axon counts") +
    theme(
        axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.3
        ),
        legend.position = "none"
    )

ggsave(plot = axon_count_level2, file.path("results", "histo", "level2_axon_counts_normal.pdf"), width = 2.2, height = 3)

# correlating normal axons with ephysio ----
cor_axon_ephysio <-
    axon_count_fascicle |>
    group_by(sample) |>
    summarize(
        axon_count = mean(normal_myelin),
        ncv_tibial_motoric = mean(ncv_tibial_motoric)
    ) |>
    mutate(log_axon_normal = log(axon_count)) |>
    ggplot(aes(x = ncv_tibial_motoric, y = axon_count)) +
    # ggplot(aes(x = ncv_tibial_motoric, y = log_axon_normal)) +
    geom_smooth(method = "lm") +
    geom_point() +
    theme_bw() +
    xlab("NCV tibial motoric (m/s)") +
    ylab("Normal axon count")

ggsave(plot = cor_axon_ephysio, file.path("results", "histo", "cor_axon_ephysio.pdf"), width = 3, height = 3)
# ggsave(plot = cor_axon_ephysio, file.path("results", "histo", "cor_log_axon_ephysio.pdf"), width = 3, height = 3)


# correlating axon counts with mySC ----
abundance_axon <-
  table(sc_merge$cluster, sc_merge$sample) |>
  as.data.frame.matrix() |>
  rownames_to_column("cell") |>
  mutate(across(where(is.numeric), function(x) x / sum(x) * 100)) |>
  pivot_longer(!cell, names_to = "sample", values_to = "count") |>
  left_join(axon_count_mean, join_by(sample)) 

abundance_axon_mySC <-
    abundance_axon |>
    dplyr::filter(cell == "mySC") |>
    ggplot(aes(x = log_axon_normal, y = count)) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_bw() +
    xlab("Log normal axon counts") +
    ylab("mySC (%)")

ggsave(file.path("results", "histo", "cor_axon_mySC.pdf"), plot = abundance_axon_mySC, width = 3, height = 3)

# abundance_axon_nmSC <-
#   abundance_axon |>
#   dplyr::filter(cell == "nmSC") |>
#   ggplot(aes(x = log_axon_normal, y = count)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   theme_bw()

# ggsave(file.path("results", "abundance", "abundance_axon_nmSC.pdf"), plot = abundance_axon_nmSC, width = 3, height = 3)

abundance_axon_repairSC <-
    abundance_axon |>
    dplyr::filter(cell == "repairSC") |>
    ggplot(aes(x = log_axon_normal, y = count)) +
    #   ggplot(aes(x = axon_normal, y = count)) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_bw() +
    xlab("Log normal axon counts") +
    ylab("repairSC (%)")

ggsave(file.path("results", "histo", "cor_axon_repairSC.pdf"), plot = abundance_axon_repairSC, width = 3, height = 3)

# abundance_axon_damageSC <-
#   abundance_axon |>
#   dplyr::filter(cell == "damageSC") |>
#   ggplot(aes(x = log_axon_normal, y = count)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   theme_bw()

# ggsave(file.path("results", "abundance", "abundance_axon_damageSC.pdf"), plot = abundance_axon_damageSC, width = 3, height = 3)

# # correlating g ratio with mySC ----
# g_ratio_mean <-
#   g_ratio |>
#   group_by(sample) |>
#   summarize(g_ratio = mean(g_ratio)) 

# abundance_g_ratio <-
#   table(sc_merge$cluster, sc_merge$sample) |>
#   as.data.frame.matrix() |>
#   rownames_to_column("cell") |>
#   mutate(across(where(is.numeric), function(x) x / sum(x) * 100)) |>
#   pivot_longer(!cell, names_to = "sample", values_to = "count") |>
#   left_join(g_ratio_mean, join_by(sample))

# abundance_g_ratio_mySC <-
#   abundance_g_ratio |>
#   dplyr::filter(cell == "mySC") |>
#   ggplot(aes(x = g_ratio, y = count)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   theme_bw()

# ggsave(file.path("results", "abundance", "abundance_g_ratio_mySC.pdf"), plot = abundance_g_ratio_mySC, width = 3, height = 3)

# abundance_g_ratio_repairSC <-
#   abundance_g_ratio |>
#   dplyr::filter(cell == "repairSC") |>
#   ggplot(aes(x = g_ratio, y = count)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   theme_bw()

# ggsave(file.path("results", "abundance", "abundance_g_ratio_repairSC.pdf"), plot = abundance_g_ratio_repairSC, width = 3, height = 3)

# axonEphysioPlot <- function(x_axis, y_axis) {
#     axon_count_fascicle |>
#         group_by(sample) |>
#         summarize(
#             x_axis = median(.data[[x_axis]]),
#             y_axis = median(.data[[y_axis]])
#         ) |>
#         ggplot(aes(x = x_axis, y = y_axis)) +
#         geom_smooth(method = "lm") +
#         geom_jitter(height = 0, width = 0.1) +
#         guides(fill = guide_legend(title = NULL)) +
#         theme_bw() +
#         xlab(x_axis) +
#         ylab(y_axis) 
# }


# x_axis_var <- axon_count_fascicle |>
#     select(cmap_ulnar:ncv_sural) |>
#     names()

# # normal axons
# normal_axon_ephysio_plots <- lapply(
#     x_axis_var,
#     function(x) {
#         axonEphysioPlot(y_axis = "normal_myelin", x_axis = x)
#     }
# )

# normal_axon_ephysio_plots_patch <- patchwork::wrap_plots(normal_axon_ephysio_plots, ncol = 3)
# ggsave(plot = normal_axon_ephysio_plots_patch, file.path("results", "histo", "normal_axon_ephysio.pdf"), width = 15, height = 15)

# # thinly myelinated axons
# thin_axon_ephysio_plots <- lapply(
#     x_axis_var,
#     function(x) {
#         axonEphysioPlot(y_axis = "thin_myelin", x_axis = x)
#     }
# )

# thin_axon_ephysio_plots_patch <- patchwork::wrap_plots(thin_axon_ephysio_plots, ncol = 3)
# ggsave(plot = thin_axon_ephysio_plots_patch, file.path("results", "histo", "thin_axon_ephysio.pdf"), width = 15, height = 15)


# # normal myelin density
# normal_density_axon_ephysio_plots <- lapply(
#     x_axis_var,
#     function(x) {
#         axonEphysioPlot(y_axis = "normal_myelin_density", x_axis = x)
#     }
# )

# normal_density_axon_ephysio_plots_patch <- patchwork::wrap_plots(normal_density_axon_ephysio_plots, ncol = 3)
# ggsave(plot = normal_density_axon_ephysio_plots_patch, file.path("results", "histo", "normal_density_axon_ephysio.pdf"), width = 15, height = 15)

# axon_plots_fascicle_ncv <- lapply(
#     c(
#         "normal_myelin",
#         "thin_myelin",
#         "very_small_axons",
#         "total_myelinated_axons",
#         "normal_myelin_density",
#         "thin_myelin_density",
#         "very_small_axons_density",
#         "total_myelinated_axons_density"
#     ),
#     function(x) {
#         axonPlotNCV(df = axon_count_fascicle, var = x)
#     }
# )

# axon_plots_fascicle_ncv_patch <- patchwork::wrap_plots(axon_plots_fascicle_ncv, ncol = 2)
# ggsave(plot = axon_plots_fascicle_ncv_patch, file.path("results", "histo", "axon_ncv_fascicle.pdf"), width = 15, height = 15)


# axon_plots_sum_ncv <- lapply(
#     c(
#         "normal_myelin",
#         "thin_myelin",
#         "very_small_axons",
#         "total_myelinated_axons",
#         "normal_myelin_density",
#         "thin_myelin_density",
#         "very_small_axons_density",
#         "total_myelinated_axons_density"
#     ),
#     function(x) {
#         axonPlotNCV(df = axon_count_sum, var = x)
#     }
# )

# axon_plots_sum_ncv_patch <- patchwork::wrap_plots(axon_plots_sum_ncv, ncol = 2)
# ggsave(plot = axon_plots_sum_ncv_patch, file.path("results", "histo", "axon_ncv_sum.pdf"), width = 15, height = 15)


# axonPlot <- function(df, y_value) {
#     df |>
#         mutate(order = as.numeric(level2)) |>
#         mutate(sample = fct_reorder(sample, order)) |>
#         ggplot(aes(x = sample, y = .data[[y_value]], fill = level2, color = level2)) +
#         geom_boxplot() +
#         geom_jitter(height = 0, width = 0.1 ) +
#         # scale_fill_manual(values = sc_merge@misc$level2_cols) +
#         scale_color_manual(values = sc_merge@misc$level2_cols) +
#         guides(fill = guide_legend(title = NULL)) +
#         theme_classic() +
#         ylab("") +
#         xlab("") +
#         ggtitle(y_value) +
#         theme(axis.text.x = element_text(
#             angle = 90,
#             hjust = 1, vjust = 0.3
#         ))
# }

# axon_plots_fascicle <- lapply(
#     c(
#         "normal_myelin",
#         "thin_myelin",
#         "very_small_axons",
#         "total_myelinated_axons",
#         "normal_myelin_density",
#         "thin_myelin_density",
#         "very_small_axons_density",
#         "total_myelinated_axons_density"
#     ),
#     function(x) {
#         axonPlot(df = axon_count_fascicle, y_value = x)
#     }
# )

# axon_plots_fascicle_patch <- patchwork::wrap_plots(axon_plots_fascicle, ncol = 2)
# ggsave(plot = axon_plots_fascicle_patch, file.path("results", "histo", "axon_count_fascicle.pdf"), width = 15, height = 15)

# axon_plots_sum <- lapply(
#     c(
#         "normal_myelin",
#         "thin_myelin",
#         "very_small_axons",
#         "total_myelinated_axons",
#         "normal_myelin_density",
#         "thin_myelin_density",
#         "very_small_axons_density",
#         "total_myelinated_axons_density"
#     ),
#     function(x) {
#         axonPlot(df = axon_count_sum, y_value = x)
#     }
# )

# axon_plots_sum_patch <- patchwork::wrap_plots(axon_plots_sum, ncol = 2)
# ggsave(plot = axon_plots_sum_patch, file.path("results", "histo", "axon_count_sum.pdf"), width = 15, height = 15)
