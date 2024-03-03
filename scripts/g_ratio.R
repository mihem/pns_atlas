# analyse g-ratio data

# load libraries ----
library(tidyverse)
library(conflicted)
library(readxl)
library(pals)
library(patchwork)
library(broom)

# read data ----
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

g_ratio <-
    read_excel(file.path("lookup", "g_ratio.xlsx")) |>
    left_join(sample_lookup, join_by(sample)) |>
    mutate(level2 = factor(level2, levels = sc_merge@misc$level2_order))

g_ratio_plot <-
    g_ratio |>
    ggplot(aes(x = sample, y = g_ratio)) +
    geom_boxplot() +
    geom_jitter(height = 0, width = 0.1, alpha = 0.3, size = 0.1) +
    theme_bw()

ggsave(file.path("results", "gratio", "g_ratio.pdf"),
    width = 10, height = 3,
    plot = g_ratio_plot
)

gratioPlot <- function(name) {
    g_ratio |>
        dplyr::filter(level2 == {{ name }}) |>
        ggplot(aes(x = axon_diameter, y = g_ratio, color = sample)) +
        geom_point(alpha = 0.5, size = 0.1) +
        geom_smooth(method = "lm", se = FALSE) +
        scale_color_manual(values = pals::cols25()) +
        theme_bw() +
        ggtitle(name)
}
g_ratio_axon_plots <- lapply(levels(g_ratio$level2), FUN = gratioPlot)
g_ratio_axon_patch <- patchwork::wrap_plots(g_ratio_level2_plots, ncol = 2)

ggsave(file.path("results", "gratio", "axon_gratio.pdf"), plot = g_ratio_axon_patch, width = 10, height = 10)

g_ratio_level2 <-
    g_ratio |>
    group_by(sample) |>
    mutate(g_ratio = mean(g_ratio)) |>
    select(sample, level2, g_ratio) |>
    distinct() |>
    ggplot(aes(x = level2, y = g_ratio, fill = level2)) +
    geom_boxplot() +
    geom_jitter(height = 0, width = 0.1) +
    scale_fill_manual(values = sc_merge@misc$level2_cols) +
    theme_classic() +
    ylab("") +
    xlab("") +
    ggtitle("g-ratio") +
    theme(
        axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.3
        ),
        legend.position = "none"
    )

ggsave(plot = g_ratio_level2, file.path("results", "gratio", "level2_g_ratio.pdf"), width = 5, height = 5)

library(dplyr)
library(purrr)

# Assuming g_ratio and axon_diameter are columns in your data frame, and it's grouped by 'sample'

lmGratio <- function(data) {
    model <- lm(g_ratio ~ axon_diameter, data = data) 
    model <- broom::tidy(model)
    result <- model$estimate[2]
    return(result)
}

models_gratio <-
    g_ratio |>
    group_by(sample) |>
    nest() |>
    mutate(slope = map_dbl(data, function(x) lmGratio(data = x))) |>
    unnest(cols = c(data)) |>
    select(sample, level2, slope) |>
    distinct()

s22_gratio <-
    g_ratio |>
    dplyr::filter(sample == "S22")

s22_gratio |>
    lm(g_ratio ~ axon_diameter, data = _) |>
    broom::tidy()


slope_level2 <-
    models_gratio |>
    ggplot(aes(x = level2, y = slope, fill = level2)) +
    geom_boxplot() +
    geom_jitter(height = 0, width = 0.1) +
    scale_fill_manual(values = sc_merge@misc$level2_cols) +
    theme_classic() +
    ylab("") +
    xlab("") +
    ggtitle("slope") +
    theme(
        axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.3
        ),
        legend.position = "none"
    )

ggsave(plot = slope_level2, file.path("results", "gratio", "level2_slope.pdf"), width = 5, height = 5)

# comparison with ncv ---
g_ratio <-
    g_ratio |>
    group_by(sample) |>
    mutate(g_ratio = mean(g_ratio),
            axon_diameter = mean(axon_diameter)) |>
    ungroup() |>
    distinct() |>
    left_join(select(models_gratio, sample, slope), join_by(sample)) |>
    distinct()

gratioEphysioPlot <- function(x_axis, y_axis) {
    g_ratio |>
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

# correlation of ephysio with gratio ---
x_axis_var <- g_ratio |>
    select(cmap_ulnar:ncv_sural) |>
    names()

gratio_ephysio_plots <- lapply(
    x_axis_var,
    function(x) {
        gratioEphysioPlot(y_axis = "g_ratio", x_axis = x)
    }
)

gratio_ephysio_plots_patch <- patchwork::wrap_plots(gratio_ephysio_plots, ncol = 3)
ggsave(plot = gratio_ephysio_plots_patch, file.path("results", "gratio", "gratio_ephysio.pdf"), width = 15, height = 15)

# correlation of ephysio with g ratio slope
gratio_slope_plots <- lapply(
    x_axis_var,
    function(x) {
        gratioEphysioPlot(y_axis = "slope", x_axis = x)
    }
)

gratio_slope_plots_patch <- patchwork::wrap_plots(gratio_slope_plots, ncol = 3)
ggsave(plot = gratio_slope_plots_patch, file.path("results", "gratio", "slope_ephysio.pdf"), width = 15, height = 15)

# correlation of ephysio with axon diameter
axon_diameter_plots <- lapply(
    x_axis_var,
    function(x) {
        gratioEphysioPlot(y_axis = "axon_diameter", x_axis = x)
    }
)

axon_diameter_plots_patch <- patchwork::wrap_plots(axon_diameter_plots, ncol = 3)
ggsave(plot = axon_diameter_plots_patch, file.path("results", "gratio", "axon_diameter_ephysio.pdf"), width = 15, height = 15)
