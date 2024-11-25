#===============================================================================
# Demographics Analysis Script 
#===============================================================================
# Purpose: Analyze and visualize patient characteristics across disease groups:
# - Age distribution
# - Nerve conduction velocity (NCV) 
# - INCAT disability scores
# - Sex distribution
#===============================================================================

# libraries ---
library(tidyverse)
library(pals)

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)

sample_lookup <-
    read_csv(file.path("lookup", "sample_lookup.csv")) |>
    janitor::clean_names() |>
    dplyr::rename(
        ncv_tibial_motoric = ncv_tibial_motoric_in_m_s,
        cmap_tibial_motoric = cmap_tibial_in_m_v,
        f_latency_tibial = min_f_latency_tibial_in_ms,
        ncv_peroneal_motoric = ncv_peroneal_motoric_in_m_s,
        cmap_peroneal_motoric = cmap_peroneal_motoric_in_m_v,
        ncv_ulnar_motoric = ncv_ulnar_motoric_in_m_s,
        cmap_ulnar = cmap_ulnar_in_m_v,
        f_latency_ulnar = min_f_latency_ulnar_in_ms,
        snap_sural = snap_sural_in_m_v,
        ncv_sural = ncv_sural_in_m_s
    ) |>
    mutate(age_calc = lubridate::time_length(difftime(nerve_date, birth_date), "years")) |>
    mutate(age_calc = floor(age_calc)) |>
    mutate(age = coalesce(age_calc, age)) |>
    dplyr::select(-age_calc) |>
    mutate(level0 = if_else(level1 == "CTRL", "CTRL", "PNP")) |>
    mutate(across(cmap_ulnar:ncv_sural, as.numeric)) |>
    mutate(level2 = factor(level2, levels = sc_merge@misc$level2_order)) |>
    mutate(incat = as.numeric(incat))

# Age Analysis ----
# Compare age distribution between disease groups using boxplots
# Statistical testing performed with scMisc::compStat()
age_plot <-
    sample_lookup |>
    ggplot(aes(x = level2, y = age, fill = level2)) +
    geom_boxplot() +
    geom_point() +
    ggsignif::geom_signif(
        comparisons = age_stats$comparisons,
        annotation = age_stats$annotation,
        textsize = 5,
        step_increase = 0.05,
        vjust = 0.7
    ) +
    theme_bw() +
    scale_fill_manual(values = sc_merge@misc$level2_cols) +
    xlab("") +
    ylab("") +
    ggtitle("age") +
    theme(legend.position = "none")

age_stats <- scMisc:::compStat(x_var = "age", group = "level2", data = sample_lookup, paired = FALSE)

ggsave(file.path("results", "demographics", "boxplot_age.pdf"), plot = age_plot, width = 5, height = 5)

# NCV Analysis   ---
# Compare tibial nerve motor conduction velocity between groups
ncv_tibial_motoric_plot <-
    sample_lookup |>
    dplyr::filter(!is.na(ncv_tibial_motoric)) |>
    ggplot(aes(x = level2, y = ncv_tibial_motoric, fill = level2)) +
    geom_boxplot() +
    geom_point() +
    ggsignif::geom_signif(
        comparisons = ncv_tibial_motoric_stats$comparisons,
        annotation = ncv_tibial_motoric_stats$annotation,
        textsize = 5,
        step_increase = 0.05,
        vjust = 0.7
    ) +
    theme_bw() +
    scale_fill_manual(values = sc_merge@misc$level2_cols) +
    xlab("") +
    ylab("") +
    ggtitle("motoric NCV tibial nerve") +
    theme(legend.position = "none")

ncv_tibial_motoric_stats <- scMisc:::compStat(x_var = "ncv_tibial_motoric", group = "level2", data = sample_lookup, paired = FALSE)

ggsave(file.path("results", "demographics", "boxplot_ncv_tibial_motoric.pdf"), plot = ncv_tibial_motoric_plot, width = 5, height = 5)

# INCAT Score Analysis ---
# Compare disability scores between disease groups 
# Scale fixed from 0-6 for consistent visualization
incat_plot <-
    sample_lookup |>
    dplyr::filter(!is.na(incat)) |>
    ggplot(aes(x = level2, y = incat, fill = level2)) +
    geom_boxplot() +
    geom_jitter(height = 0, width = 0.3) +
    ggsignif::geom_signif(
        comparisons = incat_stats$comparisons,
        annotation = incat_stats$annotation,
        textsize = 5,
        step_increase = 0.05,
        vjust = 0.7
    ) +
    theme_bw() +
    scale_fill_manual(values = sc_merge@misc$level2_cols) +
    xlab("") +
    ylab("") +
    ggtitle("INCAT score") +
    theme(legend.position = "none") + 
    ylim(0, 6)

incat_stats <- scMisc:::compStat(x_var = "incat", group = "level2", data = sample_lookup, paired = FALSE)

ggsave(file.path("results", "demographics", "boxplot_incat.pdf"), plot = incat_plot, width = 5, height = 5)

# Sex Distribution ----
# Visualize male/female distribution across disease groups
sex_plot <- 
    sample_lookup |>
    ggplot(aes(x = level2, fill = sex)) +
    geom_bar() + 
    theme_bw() +
    scale_fill_manual(values = pals::cols25(2)) +
    xlab("") +
    ylab("") +
    ggtitle("sex")

ggsave(file.path("results", "demographics", "boxplot_sex.pdf"), plot = sex_plot, width = 5, height = 5)