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
    # covs = c("sex_numeric", "age", "center_numeric") 
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
        title = "PNP", color = "Correlation"
    )

ggsave(plot = cor_pnp_plot, file.path("results", "cna", "cna_pnp_sex_age1.png"), width = 10, height = 10)

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
        title = "INCAT", color = "Correlation"
    )

ggsave(plot = cor_incat_plot, file.path("results", "cna", "cna_incat.png"), width = 10, height = 10)
ggsave(plot = cor_incat_plot, file.path("results", "cna", "cna_incat_sex_age1.png"), width = 10, height = 10)
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

    # mutate(sick = very_small_axons + thin_myelin) |>
    # mutate(ratio_normal_thin = normal_myelin / thin_myelin) |>
   # mutate(ratio_normal_sick = normal_myelin / (very_small_axons + thin_myelin)) |>

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

axon_count_median <-
    axon_count_fascicle |>
    group_by(sample) |>
    dplyr::summarize(
        axon_normal = median(normal_myelin),
    )

axon_count_median |>
dplyr::arrange(desc(axon_normal))

sc_merge@meta.data <-
    sc_merge@meta.data |>
    tibble::rownames_to_column("barcode") |>
    dplyr::left_join(axon_count_median, by = "sample") |>
    tibble::column_to_rownames(var = "barcode")

str(sc_merge@meta.data)

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

ggsave(plot = cor_axon_normal_plot, file.path("results", "cna", "cna_normal_axon1.png"), width = 10, height = 10)
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

dplyr::count(sc_merge@meta.data, ncv_tibial_motoric)

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
        title = "g ratio", color = "Correlation"
    )

ggsave(plot = cor_gratio_plot, file.path("results", "cna", "cna_gratio_sex_age.png"), width = 10, height = 10)

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
# ggsave(plot = cor_axon_diameter_plot, file.path("results", "cna", "cna_axon_diameter.png"), width = 10, height = 10)

