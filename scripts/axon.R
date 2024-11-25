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
