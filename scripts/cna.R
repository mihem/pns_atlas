# libraries ---
library(tidyverse)
library(ggthemes)
library(patchwork)
library(Seurat)
library(rcna)
library(glue)
library(RColorBrewer)
library(scales)
library(readxl)

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)
ic <- qs::qread(file.path("objects", "ic.qs"), nthread = 4)

# general settings  ----
conflicts_prefer(base::setdiff)
my_cols_25 <- pals::cols25()
my_cols_50 <- unname(Polychrome::createPalette(50, pals::cols25()))

# correlation with disease ---
sc_merge$pnp <- as.numeric(sc_merge$level2 != "CTRL")
sc_merge$sex_numeric <- as.numeric(sc_merge$sex == "male")
sc_merge$center_numeric <- as.numeric(factor(sc_merge$center))


obj <- association.Seurat(
    seurat_object = sc_merge, 
    test_var = 'pnp', 
    # test_var = 'incat_numeric', 
    samplem_key = 'sample', 
    graph_use = 'RNA_nn', 
    verbose = TRUE,
    batches = NULL, ## no batch variables to include, only works with matched design https://github.com/immunogenomics/cna/issues/11
    # covs = NULL ## no covariates to include  
    covs = c("sex_numeric", "age")
)

cor_pnp_plot <-
    FeaturePlot(
        obj,
        reduction = "umap.scvi.full",
        features = c("cna_ncorrs"),
        # features = c("cna_ncorrs_fdr10"),
        pt.size = 0.1,
        order = FALSE,
        coord.fixed = TRUE,
        raster = FALSE,
        alpha = 0.2
    ) +
    # viridis::scale_color_viridis(option = "inferno") +
    # scale_color_gradient2_tableau() +
    # scale_color_distiller(type = "div", palette = "RdBu", direction = -1) +
    scale_colour_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 0,
    ) +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA)
    ) +
    labs(
        title = "PNP", color = "Correlation", x = "UMAP1", y= "UMAP2"
    )

ggsave(plot = cor_pnp_plot, file.path("results", "cna", "cna_pnp_sex_age.png"), width = 10, height = 10)
ggsave(plot = cor_pnp_plot, file.path("results", "cna", "cna_pnp_sex_age.pdf"), width = 10, height = 10)

# immune cells
ic$pnp <- as.numeric(ic$level2 != "CTRL")
ic$sex_numeric <- as.numeric(ic$sex == "male")
ic$center_numeric <- as.numeric(factor(ic$center))

obj <- association.Seurat(
    seurat_object = ic, 
    test_var = 'pnp', 
    samplem_key = 'sample', 
    graph_use = 'RNA_nn', 
    verbose = TRUE,
    batches = NULL, ## no batch variables to include, only works with matched design https://github.com/immunogenomics/cna/issues/11
    covs = NULL
    # covs = c("sex_numeric", "age")
)

cor_pnp_plot <-
    FeaturePlot(
        obj,
        reduction = "umap.rpca",
        features = c("cna_ncorrs"),
        pt.size = 0.1,
        order = FALSE,
        coord.fixed = TRUE,
        raster = FALSE,
        alpha = 1
    ) +
    scale_colour_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 0,
    ) +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA)
    ) +
    labs(
        title = "PNP", color = "Correlation"
    )

ggsave(plot = cor_pnp_plot, file.path("results", "cna", "cna_pnp_ic_sex_age.png"), width = 10, height = 10)
ggsave(plot = cor_pnp_plot, file.path("results", "cna", "cna_pnp_ic.png"), width = 10, height = 10)

# correlation with incat ----
sc_merge$pnp <- as.numeric(sc_merge$level2 != "CTRL")

sc_merge_pnp <- subset(sc_merge, subset = level2 %in% c("CTRL"), invert = TRUE)
sc_merge_pnp$incat_numeric <- as.numeric(sc_merge_pnp$incat)
sc_merge$sex_numeric <- as.numeric(sc_merge$sex == "male")
sc_merge$center_numeric <- as.numeric(factor(sc_merge$center))

obj <- association.Seurat(
    seurat_object = sc_merge_pnp, 
    test_var = "incat_numeric", 
    samplem_key = 'sample', 
    graph_use = 'RNA_nn', 
    verbose = TRUE,
    batches = NULL , ## no batch variables to include
    covs = c("sex_numeric", "age") ## no covariates to include 
)

cor_incat_plot <-
    FeaturePlot(obj,
        reduction = "umap.scvi.full",
        features = c("cna_ncorrs"),
        pt.size = 0.1,
        order = FALSE,
        coord.fixed = TRUE,
        raster = FALSE,
        alpha = 0.2
    ) +
    scale_colour_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 0,
    ) +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA)
    ) +
    labs(
        title = "INCAT", color = "Correlation", x = "UMAP1", y= "UMAP2"
    ) 

ggsave(plot = cor_incat_plot, file.path("results", "cna", "cna_incat.png"), width = 10, height = 10)
ggsave(plot = cor_incat_plot, file.path("results", "cna", "cna_incat_sex_age.pdf"), width = 10, height = 10)
ggsave(plot = cor_incat_plot, file.path("results", "cna", "cna_incat_sex_age.png"), width = 10, height = 10)
ggsave(plot = cor_incat_plot, file.path("results", "cna", "cna_incat_sex_age_center.png"), width = 10, height = 10)

# correlation with axon count ----
# ratio to small problematic because of zeros
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

# correlation with axon counting ---
sc_merge$pnp <- as.numeric(sc_merge$level2 != "CTRL")
sc_merge$sex_numeric <- as.numeric(sc_merge$sex == "male")
sc_merge$center_numeric <- as.numeric(factor(sc_merge$center))

obj <- association.Seurat(
    seurat_object = sc_merge, 
    test_var = "axon_normal", 
    samplem_key = 'sample', 
    graph_use = 'RNA_nn', 
    verbose = TRUE,
    batches = NULL , ## no batch variables to include
    covs = NULL ## no covariates to include 
    # covs = c("sex_numeric", "age") 
)

cor_axon_normal_plot <-
    FeaturePlot(obj,
        reduction = "umap.scvi.full",
        features = c("cna_ncorrs"),
        pt.size = 0.1,
        order = FALSE,
        coord.fixed = TRUE,
        raster = FALSE,
        alpha = 0.2
    ) +
    scale_colour_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 0,
    ) +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA)
    ) +
    labs(
        title = "Normal axons", color = "Correlation"
    )

ggsave(plot = cor_axon_normal_plot, file.path("results", "cna", "cna_normal_axon.png"), width = 10, height = 10)
ggsave(plot = cor_axon_normal_plot, file.path("results", "cna", "cna_normal_axon_sex_age.png"), width = 10, height = 10)

# correlation with ncv ----
sample_lookup_ncv <- 
    read_csv(file.path("lookup", "sample_lookup.csv")) |>
    janitor::clean_names() |>
    select(sample, ncv_tibial_motoric_in_m_s)  |>
    rename(ncv_tibial_motoric = ncv_tibial_motoric_in_m_s) |>
    mutate(ncv_tibial_motoric = as.numeric(ncv_tibial_motoric))

sc_merge@meta.data <-
    sc_merge@meta.data |>
    tibble::rownames_to_column("barcode") |>
    dplyr::left_join(sample_lookup_ncv, by = "sample") |>
    tibble::column_to_rownames(var = "barcode")

samples_no_ncv <-
    dplyr::count(sc_merge@meta.data, ncv_tibial_motoric, sample) |>
    dplyr::filter(is.na(ncv_tibial_motoric)) |>
    pull(sample)

sc_merge_ncv <- subset(sc_merge, subset = sample %in% samples_no_ncv, invert = TRUE)

obj <- association.Seurat(
    seurat_object = sc_merge_ncv, 
    test_var = "ncv_tibial_motoric", 
    samplem_key = 'sample', 
    graph_use = 'RNA_nn', 
    verbose = TRUE,
    batches = NULL , ## no batch variables to include
    covs = c("sex_numeric", "age", "center_numeric") ## no covariates to include 
)

cor_ncv_plot <-
    FeaturePlot(obj,
        reduction = "umap.scvi.full",
        features = c("cna_ncorrs"),
        pt.size = 0.1,
        order = FALSE,
        coord.fixed = TRUE,
        raster = FALSE,
        alpha = 0.2
    ) +
    scale_colour_gradient2(
        low = muted("blue"),
        mid = "white",
        high = muted("red"),
        midpoint = 0,
    ) +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA)
    ) +
    labs(
        title = "NCV tibial motoric", color = "Correlation"
    )

ggsave(plot = cor_ncv_plot, file.path("results", "cna", "cna_ncv_tibial_motoric_cov.png"), width = 10, height = 10)

# correlation with g ratio ----
sc_merge@meta.data <-
    sc_merge@meta.data |>
    tibble::rownames_to_column("barcode") |>
    dplyr::left_join(select(g_ratio, sample, g_ratio, axon_diameter, slope), by = "sample") |>
    rename(g_ratio_slope = slope) |>
    tibble::column_to_rownames(var = "barcode")

ic@meta.data <-
    ic@meta.data |>
    tibble::rownames_to_column("barcode") |>
    dplyr::left_join(select(g_ratio, sample, g_ratio, axon_diameter, slope), by = "sample") |>
    rename(g_ratio_slope = slope) |>
    tibble::column_to_rownames(var = "barcode")

obj <- association.Seurat(
    seurat_object = sc_merge, 
    test_var = "g_ratio", 
    samplem_key = 'sample', 
    graph_use = 'RNA_nn', 
    verbose = TRUE,
    batches = NULL, ## no batch variables to include
    covs = c("sex_numeric", "age") ## no covariates to include 
)

cor_gratio_plot <-
    FeaturePlot(
        obj,
        reduction = "umap.scvi.full",
        features = c("cna_ncorrs"),
        # features = c("cna_ncorrs_fdr10"),
        pt.size = 0.1,
        order = FALSE,
        coord.fixed = TRUE,
        raster = FALSE,
        alpha = 0.2
    ) +
    scale_colour_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 0,
    ) +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA)
    ) +
    labs(
        title = "g ratio", color = "Correlation", x = "UMAP1", y = "UMAP2"
    )

ggsave(plot = cor_gratio_plot, file.path("results", "cna", "cna_gratio_sex_age.png"), width = 10, height = 10)
ggsave(plot = cor_gratio_plot, file.path("results", "cna", "cna_gratio_sex_age.pdf"), width = 10, height = 10)

# correlation with axon diameter ----
obj <- association.Seurat(
    seurat_object = sc_merge, 
    test_var = "axon_diameter", 
    samplem_key = 'sample', 
    graph_use = 'RNA_nn', 
    verbose = TRUE,
    batches = NULL, ## no batch variables to include
    # covs = NULL
    covs = c("sex_numeric", "age") ## no covariates to include 
)

cor_axon_diameter_plot <-
    FeaturePlot(
        obj,
        reduction = "umap.scvi.full",
        features = c("cna_ncorrs"),
        # features = c("cna_ncorrs_fdr10"),
        pt.size = 0.1,
        order = FALSE,
        coord.fixed = TRUE,
        raster = FALSE,
        alpha = 0.2
    ) +
    scale_colour_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 0,
    ) +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA)
    ) +
    labs(
        title = "axon diameter", color = "Correlation"
    )

ggsave(plot = cor_axon_diameter_plot, file.path("results", "cna", "cna_axon_diameter_sex_age.png"), width = 10, height = 10)

