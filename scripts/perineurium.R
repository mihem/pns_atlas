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

# meta data 
perineurium_lookup <- read_csv(file.path("lookup", "perineurium_lookup.csv"))

perineurium <-
    read_csv(file.path("lookup", "perineurium_measurement.csv")) |>
    dplyr::filter(is.na(remove) | !remove == "yes") |>
    left_join(perineurium_lookup) |>
    mutate(diagnosis = factor(diagnosis, levels = diagnosis_order)) |>
    mutate(sample_ordered = reorder(sample, as.numeric(diagnosis)))  |>
    dplyr::mutate(outer_inner_diameter_ratio = outer_diameter / inner_diameter) |>
    dplyr::mutate(outer_inner_area_ratio = outer_area / inner_area) 

# sanity check
all(perineurium$inner_area < perineurium$outer_area)
all(perineurium$inner_diameter < perineurium$outer_diameter)
all(!duplicated(count(perineurium, sample, diagnosis)$sample))
perineurium[perineurium$inner_diameter > perineurium$outer_diameter, ]
print(dplyr::count(perineurium, sample, diagnosis), n = Inf)
select(arrange(perineurium, desc(outer_inner_area_ratio)), sample, outer_area, inner_area, outer_inner_area_ratio)


# plot perineurium diameter ratio ----
perineurium_outer_inner_diameter <-
    perineurium |>
    ggplot(aes(x = sample_ordered, y = outer_inner_diameter_ratio, fill = diagnosis)) +
    geom_boxplot() +
    geom_point() +
    theme_classic() +
    scale_fill_manual(values = diagnosis_col) 

ggsave(
    file.path("results", "perineurium", "perineurium_outer_inner_diameter_ratio.pdf"),
    plot = perineurium_outer_inner_diameter,
    width = 13,
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
    width = 13,
    height = 5
)

perineurium_outer_inner_area_ratio <-
    perineurium |>
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
    width = 12,
    height = 5
)
