# ===============================================================================
# Analysis of perineurium
# ===============================================================================
# Purpose: Analyze manual measurements of perineurium thickness
# ===============================================================================

# load libraries ----
library(tidyverse)
library(readxl)

# read data ----
diagnosis_order <- c("CTRL", "VN", "CIDP", "CIAP", "PPN", "DPN", "OIN", "ONIN", "MNC", "IN")
diagnosis_col <- setNames(pals::cols25(length(diagnosis_order)), diagnosis_order)

sample_lookup <-
    read_csv(file.path("lookup", "sample_lookup.csv")) |>
    dplyr::select(sample, level1, level2, center, INCAT)

perineurium <-
    read_csv(file.path("lookup", "perineurium_measurement.csv")) |>
    dplyr::filter(remove == "no") |>
    mutate(original_name = sub("_EMA\\.svs \\[0\\]", "", x = original_name)) |>
    left_join(sample_lookup, join_by(sample)) |>
    mutate(center = coalesce(center, center_CH)) |>
    mutate(diagnosis = coalesce(level2, diagnosis_CH)) |>
    mutate(diagnosis = factor(diagnosis, levels = diagnosis_order)) |>
    mutate(sample_ordered = reorder(sample, as.numeric(diagnosis))) |>
    mutate(area = outer_area - inner_area)

# sanity check
all(perineurium$inner_area < perineurium$outer_area)
dplyr::count(perineurium, diagnosis, sample)
unique(perineurium$sample)
dplyr::count(perineurium, sample, original_name)

# plot perineurium diameter ratio ----
perineurium_outer_inner_diameter <-
    perineurium |>
    dplyr::filter(!is.na(outer_diameter)) |>
    dplyr::mutate(outer_inner_diameter_ratio = outer_diameter / inner_diameter) |>
    ggplot(aes(x = sample, y = outer_inner_diameter_ratio, fill = level2)) +
    geom_boxplot() +
    geom_point() +
    theme_classic() +
    scale_fill_manual(values = diagnosis_col) 

ggsave(
    file.path("results", "perineurium", "perineurium_outer_inner_diameter_ratio.pdf"),
    plot = perineurium_outer_inner_diameter,
    width = 5,
    height = 5
)

perineurium_inner_outer_diameter <-
    perineurium |>
    dplyr::filter(!is.na(outer_diameter)) |>
    dplyr::mutate(inner_outer_diameter_ratio = inner_diameter / outer_diameter) |>
    ggplot(aes(x = sample, y = inner_outer_diameter_ratio, fill = level2)) +
    geom_boxplot() +
    geom_point() +
    theme_classic() +
    scale_fill_manual(values = diagnosis_col)

ggsave(
    file.path("results", "perineurium", "perineurium_inner_outer_diameter_ratio.pdf"),
    plot = perineurium_inner_outer_diameter,
    width = 5,
    height = 5
)

# plot perineurial area ----
perineurium_area <-
    perineurium |>
    ggplot(aes(x = sample_ordered,
               y = area, 
               fill = diagnosis)) +
    geom_boxplot() +
    geom_point() +
    theme_classic() + 
    scale_fill_manual(values = diagnosis_col)

ggsave(
    file.path("results", "perineurium", "perineurium_area.pdf"),
    plot = perineurium_area,
    width = 8,
    height = 5
)

perineurium_outer_inner_area_ratio <-
    perineurium |>
    dplyr::mutate(outer_inner_area_ratio = outer_area / inner_area) |>
    ggplot(aes(
      x = sample_ordered,
      y = outer_inner_area_ratio,
      fill = diagnosis
    )) +
    geom_boxplot() +
    geom_point() +
    theme_classic() +
    scale_fill_manual(values = diagnosis_col)

ggsave(
    file.path("results", "perineurium", "perineurium_outer_inner_area_ratio.pdf"),
    plot = perineurium_outer_inner_area_ratio,
    width = 8,
    height = 5
)
