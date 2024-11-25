#===============================================================================
# G-Ratio Analysis Script
#===============================================================================
# Purpose: Analyze and visualize g-ratio measurements from nerve biopsies and 
# correlate with electrophysiological parameters
#
# Methods:
# 1. Load and preprocess g-ratio and sample metadata
# 2. Generate boxplots of g-ratio and axon diameter distributions
# 3. Analyze correlations between g-ratio/axon diameter and nerve conduction data
#
# Key Variables:
# - g_ratio: Ratio of inner to outer diameter of myelinated nerve fibers
# - axon_diameter: Diameter of the axon (µm)
# - ncv_*: Nerve conduction velocities for different nerves
# - cmap_*: Compound muscle action potentials
#===============================================================================

# load libraries ----
library(tidyverse)
library(conflicted)
library(readxl)
library(pals)
library(patchwork)
library(broom)

# Load and preprocess sample metadata ---
# Includes clinical parameters and nerve conduction studies
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

# Load and join g-ratio measurements with sample metadata ----
g_ratio <-
    read_excel(file.path("lookup", "g_ratio.xlsx")) |>
    left_join(sample_lookup, join_by(sample)) |>
    mutate(level2 = factor(level2, levels = sc_merge@misc$level2_order))

# Generate boxplots showing g-ratio distribution by sample ----
g_ratio_plot <-
    g_ratio |>
    ggplot(aes(x = sample, y = g_ratio)) +
    geom_boxplot() +
    geom_jitter(height = 0, width = 0.1, alpha = 0.3, size = 0.1) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    xlab("") + 
    ylab("") +
    ggtitle("g-ratio") 

ggsave(file.path("results", "gratio", "g_ratio.pdf"),
    width = 6, height = 3,
    plot = g_ratio_plot
)

# Generate boxplots of g-ratio by disease group (level2) ----
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

ggsave(plot = g_ratio_level2, file.path("results", "gratio", "level2_g_ratio.pdf"), width = 2, height = 3)

# axon diameter level2 ----
axon_diameter_level2 <-
    g_ratio |>
    group_by(sample) |>
    mutate(axon_diameter = mean(axon_diameter)) |>
    select(sample, level2, axon_diameter) |>
    distinct() |>
    ggplot(aes(x = level2, y = axon_diameter, fill = level2)) +
    geom_boxplot() +
    geom_jitter(height = 0, width = 0.1) +
    scale_fill_manual(values = sc_merge@misc$level2_cols) +
    theme_classic() +
    ylab("") +
    xlab("") +
    ggtitle("Axon diameter (µm)") +
    theme(
        axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.3
        ),
        legend.position = "none"
    )
ggsave(plot = axon_diameter_level2, file.path("results", "gratio", "level2_axon_diameter.pdf"), width = 2, height = 3)

# Analyze correlation between g-ratio and nerve conduction velocity ----
cor_g_ratio_ephysio <-
    g_ratio |>
    group_by(sample) |>
    summarize(
        g_ratio = mean(g_ratio),
        ncv_tibial_motoric = mean(ncv_tibial_motoric)
    ) |>
    ggplot(aes(x = ncv_tibial_motoric, y = g_ratio)) +
    geom_smooth(method = "lm") +
    geom_point() + 
    theme_bw() +
    xlab("NCV tibial motoric (m/s)") +
    ylab("Normal axon count")
ggsave(plot = cor_g_ratio_ephysio, file.path("results", "gratio", "cor_g_ratio_ephysio.pdf"), width = 3, height = 3)

# Analyze correlation between axon diameter and nerve conduction velocity ---- 
cor_axon_diameter_ephysio <-
    g_ratio |>
    group_by(sample) |>
    summarize(
        axon_diameter = mean(axon_diameter),
        ncv_tibial_motoric = mean(ncv_tibial_motoric)
    ) |>
    ggplot(aes(x = ncv_tibial_motoric, y = axon_diameter)) +
    geom_smooth(method = "lm") +
    geom_point() + 
    theme_bw() +
    xlab("NCV tibial motoric (m/s)") +
    ylab("Axon diameter (µm)")
ggsave(plot = cor_axon_diameter_ephysio, file.path("results", "gratio", "cor_axon_diameter_ephysio.pdf"), width = 3, height = 3)

sc_merge@meta.data |>
    dplyr::select(sample, age, axon_normal, g_ratio, axon_diameter) |>
    distinct() |>
    write_csv(file.path("results", "qc", "meta_data.csv"))
