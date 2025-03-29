# ===============================================================================
# Figure Reproducibility Preparation Script
# ===============================================================================
# Purpose: Prepare Seurat object for qmd file (reproducibility the figures)
# ===============================================================================

# Load necessary libraries
library(Seurat)
library(qs)
library(tidyverse)
library(miloDE)
library(SingleCellExperiment)

source(file.path("scripts", "dotplot_functions.R"))
source(file.path("scripts", "vlnplot_functions.R"))

# Load data
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)
ic <- qs::qread(file.path("objects", "ic.qs"), nthread = 4)

# Figure 1 ----
# Prepare Seurat object for reproducibility
# DietSeurat reduces the size of the Seurat object by keeping only the necessary data
umap_figure <- DietSeurat(
    sc_merge,
    counts = TRUE,
    data = FALSE,
    scale.data = FALSE,
    assays = "RNA",
    dimreducs = c("umap.scvi.full")
)

# Remove unnecessary data to further reduce the object size
umap_figure$RNA$counts <- NULL
umap_figure$RNA$data <- NULL
umap_figure$RNA$scale.data <- Matrix::Matrix(
    0,
    nrow = nrow(umap_figure$RNA),
    ncol = ncol(umap_figure$RNA)
)

umap_figure@meta.data <-
    umap_figure@meta.data |>
    tibble::rownames_to_column("barcode") |>
    dplyr::select(
        barcode,
        sample,
        center,
        sex,
        age,
        incat,
        milbrandt_sciatic_label_full,
        milbrandt_sciatic_label_full.score,
        suter_p60_label_full,
        suter_p60_label_full.score,
        nCount_RNA,
        nFeature_RNA,
        percent_mt,
        scDblFinder.score
    ) |>
    tibble::column_to_rownames("barcode")

umap_figure@commands <- list()
umap_figure@tools <- list()

# store colors in Seurat object
set.seed(123)
my_cols_50 <- unname(Polychrome::createPalette(50, pals::cols25()))
umap_figure@misc$sample_cols <- my_cols_50[1:37]

# store markers in Seurat object
markers_dotplot <- read_csv(file.path("lookup", "markers.csv")) |>
    select(dotplot_jolien) |>
    drop_na() |>
    pull()

umap_figure@misc$markers_dotplot <- markers_dotplot

# Save the processed Seurat object to a file
qs::qsave(umap_figure, file.path("docs", "umap_figure.qs"))

# Figure 2 ----

# Remove unnecessary data to further reduce the object size
ic_figure <- DietSeurat(
    ic,
    counts = TRUE,
    data = FALSE,
    scale.data = FALSE,
    assays = "RNA",
    dimreducs = c("umap.rpca")
)

# Remove unnecessary data to further reduce the object size
ic_figure$RNA$counts <- NULL
ic_figure$RNA$scale.data <- NULL

ic_figure@meta.data <- ic_figure@meta.data["stroke_label"]
ic_figure@commands <- list()
ic_figure@tools <- list()


# store markers in Seurat object
markers_dotplot_ic <- read_csv(file.path("lookup", "markers.csv")) |>
    select(dotplot_ic_short) |>
    drop_na() |>
    pull()

ic_figure@misc$markers_dotplot_ic <- markers_dotplot_ic

qsave(ic_figure, file.path("docs", "ic_figure.qs"))

## prepare data for dotplot IC SAMC ----
samc_genes <- c("SPP1", "APOE", "LPL", "FABP5", "GPNMB")

dotplot_data_ic_samc <-
    DotPlotData(
        object = ic,
        features = samc_genes,
        dot.min = 0.01,
    )

qsave(dotplot_data_ic_samc, file.path("docs", "dotplot_data_ic_samc.qs"))

# subset the object to only include the B and Plasma clusters
b_plasma <- subset(ic, idents = c("Plasma", "B"))

# Remove unnecessary data to further reduce the object size
b_plasma_figure <- DietSeurat(
    b_plasma,
    counts = TRUE,
    data = FALSE,
    scale.data = FALSE,
    assays = "RNA",
    dimreducs = NULL
)

b_plasma_figure@meta.data <-
    b_plasma_figure@meta.data |>
    dplyr::select(ic_cluster)

b_plasma_genes <- c(
    "IGHM",
    "IGHD",
    "IGHG1",
    "IGHG2",
    "IGHG3",
    "IGHG4",
    "IGHA1",
    "IGHA2"
)

dotplot_data_b_plasma <-
    DotPlotData(
        object = b_plasma,
        features = b_plasma_genes,
        dot.min = 0.01,
        scale = FALSE
    )

qsave(dotplot_data_b_plasma, file.path("docs", "dotplot_data_b_plasma.qs"))

# save enrichment results as qs object
enrichr_macro18 <- readxl::read_excel(
    file.path("results", "enrichr", "enrichr_Macro18.xlsx"),
    sheet = "GO_Biological_Process_2023"
)
qsave(enrichr_macro18, file.path("docs", "enrichr_macro18.qs"))

# Figure 3 A ----
propeller_PNP_CTRL <-
    scMisc::propellerCalc(
        seu_obj1 = sc_merge,
        condition1 = "PNP",
        condition2 = "CTRL",
        cluster_col = "cluster",
        meta_col = "level0",
        lookup = sample_lookup,
        sample_col = "sample",
        formula = "~0 + level0",
        min_cells = 30
    ) |>
    dplyr::filter(abs(log2ratio) > 0.5)

qsave(propeller_PNP_CTRL, file.path("docs", "propeller_PNP_CTRL.qs"))

propeller_PNP_CTRL_ic <-
    scMisc::propellerCalc(
        seu_obj1 = ic,
        condition1 = "PNP",
        condition2 = "CTRL",
        cluster_col = "ic_cluster",
        meta_col = "level0",
        lookup = sample_lookup,
        sample_col = "sample",
        formula = "~0 + level0",
        min_cells = 30
    ) |>
    dplyr::filter(abs(log2ratio) > 0.5)

qsave(propeller_PNP_CTRL_ic, file.path("docs", "propeller_PNP_CTRL_ic.qs"))

# Figure 3B CNA  ----
# PNP
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)
sc_merge$pnp <- as.numeric(sc_merge$level2 != "CTRL")
sc_merge$sex_numeric <- as.numeric(sc_merge$sex == "male")
sc_merge$center_numeric <- as.numeric(factor(sc_merge$center))

cna_pnp <- association.Seurat(
    seurat_object = sc_merge,
    test_var = 'pnp',
    samplem_key = 'sample',
    graph_use = 'RNA_nn',
    verbose = TRUE,
    batches = NULL, ## no batch variables to include, only works with matched design https://github.com/immunogenomics/cna/issues/11
    covs = c("sex_numeric", "age")
)
cna_pnp$cna_ncorrs_pnp <- cna_pnp$cna_ncorrs

# g ratio
cna_gratio <- association.Seurat(
    seurat_object = sc_merge,
    test_var = "g_ratio",
    samplem_key = 'sample',
    graph_use = 'RNA_nn',
    verbose = TRUE,
    batches = NULL, ## no batch variables to include
    covs = c("sex_numeric", "age")
)

cna_pnp$cna_ncorrs_gratio <- cna_gratio$cna_ncorrs

# Remove unnecessary data to further reduce the object size
cna_pnp_gratio_figure <- DietSeurat(
    cna_pnp,
    counts = TRUE,
    data = FALSE,
    scale.data = FALSE,
    assays = "RNA",
    dimreducs = c("umap.scvi.full")
)

cna_pnp_gratio_figure$RNA$counts <- NULL
cna_pnp_gratio_figure$RNA$scale.data <- NULL
cna_pnp_gratio_figure@commands <- list()
cna_pnp_gratio_figure@tools <- list()
cna_pnp_gratio_figure@meta.data <-
    cna_pnp_gratio_figure@meta.data |>
    tibble::rownames_to_column("barcode") |>
    dplyr::select(barcode, cna_ncorrs_pnp, cna_ncorrs_gratio) |>
    tibble::column_to_rownames("barcode")

qsave(cna_pnp_gratio_figure, file.path("docs", "cna_pnp_gratio_figure.qs"))

# Prepare cna_incat
sc_merge_pnp <- subset(sc_merge, subset = level2 %in% c("CTRL"), invert = TRUE)
sc_merge_pnp$incat_numeric <- as.numeric(sc_merge_pnp$incat)
sc_merge$sex_numeric <- as.numeric(sc_merge$sex == "male")

cna_incat <- association.Seurat(
    seurat_object = sc_merge_pnp,
    test_var = 'incat_numeric',
    samplem_key = 'sample',
    graph_use = 'RNA_nn',
    verbose = TRUE,
    batches = NULL,
    covs = c("sex_numeric", "age")
)

# Remove unnecessary data to further reduce the object size
cna_incat_figure <- DietSeurat(
    cna_incat,
    counts = TRUE,
    data = FALSE,
    scale.data = FALSE,
    assays = "RNA",
    dimreducs = c("umap.scvi.full")
)

cna_incat_figure$RNA$counts <- NULL
cna_incat_figure$RNA$scale.data <- NULL
cna_incat_figure@commands <- list()
cna_incat_figure@tools <- list()
cna_incat_figure@meta.data <-
    cna_incat_figure@meta.data |>
    tibble::rownames_to_column("barcode") |>
    dplyr::select(barcode, cna_ncorrs) |>
    tibble::column_to_rownames("barcode")

qsave(cna_incat_figure, file.path("docs", "cna_incat_figure.qs"))

# Figure 3C ----
sheets <- readxl::excel_sheets(
    path = file.path("results", "de", "pnp_ctrl_pseudobulk.xlsx")
)
cl_sig <-
    lapply(
        sheets,
        function(sheet) {
            readxl::read_xlsx(
                path = file.path("results", "de", "pnp_ctrl_pseudobulk.xlsx"),
                sheet = sheet
            ) |>
                dplyr::filter(p_val_adj < 0.05) |>
                nrow()
        }
    )
result <- tibble(
    cluster = sheets,
    n = unlist(cl_sig),
)

qs::qsave(result, file.path("docs", "pnp_ctrl_pseudobulk_de.qs"))

# Figure 3D and 4G ----
milo_DE <- qs::qread(file.path("objects", "milo_DE.qs"), nthread = 4)

# Print size of each slot in the Milo object
slots <- slotNames(milo_DE)
for (slot in slots) {
    size <- object.size(slot(milo_DE, slot))
    print(sprintf("Slot %s: %s", slot, format(size, units = "auto")))
}

# Remove unnecessary data to reduce the object size
empty_matrix <- Matrix::Matrix(
    0,
    nrow = nrow(milo_DE),
    ncol = ncol(milo_DE),
    sparse = TRUE
)
rownames(empty_matrix) <- rownames(milo_DE)
colnames(empty_matrix) <- colnames(milo_DE)

for (i in assayNames(milo_DE)) {
    assay(milo_DE, i) <- empty_matrix
}

colData(milo_DE) <- NULL
milo_DE@graph <- list()
reducedDims(milo_DE) <- reducedDims(milo_DE)["UMAP.SCVI.FULL"]

# Print the size of the object
print(object.size(milo_DE), units = "Gb")

de_stat <-
    lapply(
        c("pnp", "vn", "cidp", "ciap"),
        function(condition) {
            qs::qread(file.path(
                "objects",
                paste0("milo_de_stat_", condition, "_ctrl.qs")
            ))
        }
    )

stat_de_magnitude <-
    lapply(
        de_stat,
        rank_neighbourhoods_by_DE_magnitude
    )
names(stat_de_magnitude) <- c("pnp", "vn", "cidp", "ciap")

milo_figure <- list(
    obj = milo_DE,
    stat = stat_de_magnitude
)

print(object.size(milo_figure$obj), units = "GB")
print(object.size(milo_figure$stat), units = "GB")

qsave(milo_figure, file.path("docs", "milo_figure.qs"))

# Figure 3E ----
# calculate DE PNP vs Ctrl for each cluster
dePseudo <- function(seu_obj, cell_type_col, label_col) {
    res <- Libra::run_de(
        seu_obj,
        replicate_col = "sample",
        cell_type_col = cell_type_col,
        label_col = label_col,
        min_cells = 3,
        min_reps = 2,
        min_feature = 0,
        de_family = "pseudobulk",
        de_method = "edgeR",
        de_type = "LRT",
        n_threads = 6
    )
    res <- arrange(res, cell_type, desc(avg_logFC))
    res_split <- split(res, res$cell_type)
    return(res_split)
}

de_pseudo_pnp_ctrl <- dePseudo(
    sc_merge,
    cell_type_col = "cluster",
    label_col = "level0"
)
qsave(de_pseudo_pnp_ctrl, file.path("docs", "de_pseudo_pnp_ctrl.qs"))

qsave(de_pseudo_pnp_ctrl, file.path("docs", "de_pseudo_pnp_ctrl.qs"))

# Figure 4A ----

preStackedPlot <- function(object, x_axis, y_axis, x_order, y_order) {
    result_wide <- mutate(
        rownames_to_column(
            as.data.frame.matrix(table(
                object@meta.data[[y_axis]],
                object@meta.data[[x_axis]]
            )),
            "cell"
        ),
        across(where(is.numeric), function(x) x / sum(x) * 100)
    )
    result_long <- dplyr::filter(
        mutate(
            mutate(
                pivot_longer(
                    result_wide,
                    !cell,
                    names_to = "type",
                    values_to = "count"
                ),
                cell = factor(cell, levels = y_order)
            ),
            type = factor(type, levels = x_order)
        ),
        count != 0
    )
    return(result_long)
}

abundance_main_clusters <-
    preStackedPlot(
        object = sc_merge,
        x_axis = "level2",
        y_axis = "cluster",
        x_order = sc_merge@misc$level2_order,
        y_order = sc_merge@misc$cluster_order
    )

qsave(abundance_main_clusters, file.path("docs", "abundance_main_clusters.qs"))

# Figure 4B ----

## PNP subtype main clusters
propeller_PNP_subtypes_main <-
    lapply(
        c("CIDP", "CIAP", "VN"),
        function(condition) {
            scMisc::propellerCalc(
                seu_obj1 = sc_merge,
                condition1 = condition,
                condition2 = "CTRL",
                cluster_col = "cluster",
                meta_col = "level2",
                lookup = sample_lookup,
                sample_col = "sample",
                formula = "~0 + level2",
                min_cells = 30
            ) |>
                dplyr::filter(abs(log2ratio) > 0.5)
        }
    )

names(propeller_PNP_subtypes_main) <- c("CIDP", "CIAP", "VN")
qsave(
    propeller_PNP_subtypes_main,
    file.path("docs", "propeller_PNP_subtypes_main.qs")
)

## PNP subtype immune cell clusters
propeller_PNP_subtypes_ic <-
    lapply(
        c("CIDP", "CIAP", "VN"),
        function(condition) {
            scMisc::propellerCalc(
                seu_obj1 = ic,
                condition1 = condition,
                condition2 = "CTRL",
                cluster_col = "cluster",
                meta_col = "level2",
                lookup = sample_lookup,
                sample_col = "sample",
                formula = "~0 + level2",
                min_cells = 30
            ) |>
                dplyr::filter(abs(log2ratio) > 0.5)
        }
    )

names(propeller_PNP_subtypes_main) <- c("CIDP", "CIAP", "VN")
qsave(
    propeller_PNP_subtypes_main,
    file.path("docs", "propeller_PNP_subtypes_main.qs")
)

# Figure 4F ----
quanti <-
    read_csv(file.path("lookup", "xenium_manual_quantification.csv")) |>
    pivot_longer(!sample, names_to = "name", values_to = "value") |>
    separate_wider_delim(name, delim = "_", names = c("variable", "area"))

quanti_stats <-
    quanti |>
    dplyr::filter(area != "all") |>
    group_by(sample, variable) |>
    summarize(
        mean = mean(value, na.rm = TRUE),
        sum = sum(value, na.rm = TRUE)
    ) |>
    pivot_longer(c(mean, sum), names_to = "name", values_to = "value") |>
    unite("variable", variable, name)

quanti_result <-
    quanti |>
    dplyr::filter(area == "all") |>
    select(sample, variable, value) |>
    bind_rows(quanti_stats) |>
    arrange(sample) |>
    pivot_wider(names_from = variable, values_from = value) |>
    left_join(select(sample_lookup, sample, level2), join_by(sample)) |>
    mutate(level2 = factor(level2, levels = sc_merge@misc$level2_order)) |>
    mutate(epiarea = nervearea - endoarea_sum) |>
    mutate(epiPTPRC = nervePTPRC - endoPTPRC_sum) |>
    mutate(epiMS4A1 = nerveMS4A1 - endoMS4A1_sum) |>
    mutate(epiCD3E = nerveCD3E - endoCD3E_sum) |>
    mutate(epiPTPRC_density = epiPTPRC / epiarea) |>
    mutate(epiMS4A1_density = epiMS4A1 / epiarea) |>
    mutate(epiCD3E_density = epiCD3E / epiarea) |>
    mutate(endoPTPRC_density_sum = endoPTPRC_sum / endoarea_sum) |>
    mutate(endoMS4A1_density_sum = endoMS4A1_sum / endoarea_sum) |>
    mutate(endoCD3E_density_sum = endoCD3E_sum / endoarea_sum)

qsave(quanti_result, file.path("docs", "manual_xenium_quantification.qs"))

# Figure 5
perineurial <- subset(
    sc_merge,
    cluster %in%
        c("periC1", "periC2", "periC3") &
        level2 %in% c("CTRL", "CIDP", "VN", "CIAP")
)
perineurial$level2 <- factor(
    perineurial$level2,
    levels = c("CTRL", "CIDP", "VN", "CIAP")
)

# Remove unnecessary data to further reduce the object size
perineurial_figure <- DietSeurat(
    perineurial,
    counts = FALSE,
    data = TRUE,
    scale.data = FALSE,
    assays = "RNA",
    dimreducs = NULL
)


# Remove unnecessary data to further reduce the object size
perineurial_figure$RNA$counts <- NULL
perineurial_figure$RNA$scale.data <- NULL

perineurial_figure@meta.data <- NULL
perineurial_figure@commands <- list()
perineurial_figure@tools <- list()

slots <- slotNames(perineurial_figure)
for (slot in slots) {
    size <- object.size(slot(perineurial_figure, slot))
    print(paste0("Slot ", slot, ": ", format(size, units = "auto")))
}

qsave(perineurial_figure, file.path("docs", "perineurial_figure.qs"))


# Supplementary Figure 1 ----
## Supplementary Figure 1A ---
sample_lookup <-
    read_csv(file.path("lookup", "sample_lookup.csv")) |>
    janitor::clean_names() |>
    dplyr::rename(
        ncv_tibial_motoric = ncv_tibial_motoric_in_m_s,
        cmap_tibial_motoric = cmap_tibial_in_m_v,
        f_latency_tibial = min_f_latency_tibial_in_ms,
        ncv_peroneal_motoric = ncv_peroneal_motoric_in_m_s,
        ncv_ulnar_motoric = ncv_ulnar_motoric_in_m_s,
        cmap_ulnar = cmap_ulnar_in_m_v,
        f_latency_ulnar = min_f_latency_ulnar_in_ms,
        snap_sural = snap_sural_in_m_v,
        ncv_sural = ncv_sural_in_m_s
    ) |>
    mutate(
        age_calc = lubridate::time_length(
            difftime(nerve_date, birth_date),
            "years"
        )
    ) |>
    mutate(age_calc = floor(age_calc)) |>
    mutate(age = coalesce(age_calc, age)) |>
    dplyr::select(-age_calc) |>
    mutate(level0 = if_else(level1 == "CTRL", "CTRL", "PNP")) |>
    mutate(across(cmap_ulnar:ncv_sural, as.numeric)) |>
    mutate(level2 = factor(level2, levels = sc_merge@misc$level2_order)) |>
    mutate(incat = as.numeric(incat))

demographics <-
    sample_lookup |>
    dplyr::select(sample, age, sex, level2, incat, ncv_tibial_motoric)


qs::qsave(demographics, file.path("docs", "demographics.qs"))

abundance_main_clusters_sample <-
    preStackedPlot(
        object = sc_merge,
        x_axis = "sample",
        y_axis = "cluster",
        x_order = unique(sc_merge$sample),
        y_order = sc_merge@misc$cluster_order
    )


qsave(
    abundance_main_clusters_sample,
    file.path("docs", "abundance_main_clusters_sample.qs")
)

# Supplementary Figure 2 ----
## Supplementary Figure 2A ----
# load modified Seurat data

Idents(sc_merge) <- factor(
    sc_merge$cluster,
    levels = rev(sc_merge@misc$cluster_order)
)

dotplot_data <-
    DotPlotData(
        object = sc_merge,
        features = umap_figure@misc$markers_dotplot,
        dot.min = 0.01,
    )

qsave(dotplot_data, file.path("docs", "dotplot_data.qs"))

## Supplementary Figure 2B ----
enrichr_clusters <- c("periC1", "periC2", "periC3")
names(enrichr_clusters) <- c("periC1", "periC2", "periC3")

enrichr_periC <-
    lapply(
        enrichr_clusters,
        function(condition) {
            readxl::read_excel(
                file.path(
                    "results",
                    "enrichr",
                    paste0("enrichr_", condition, ".xlsx")
                ),
                sheet = "GO_Biological_Process_2023"
            )
        }
    )

qsave(enrichr_periC, file.path("docs", "enrichr_periC.qs"))

## Supplementary Figure 2C
ec <- subset(sc_merge, subset = RNA_snn_res.0.7 %in% c("7", "10", "11", "19"))
ec_rosmap <- subset(
    ec,
    subset = rosmap_label %in% c("aEndo", "capEndo", "vEndo")
)

ec_rosmap <- DietSeurat(
    ec_rosmap,
    counts = TRUE,
    data = FALSE,
    scale.data = FALSE,
    assays = "RNA",
    dimreducs = c("umap.scvi.full")
)

# Remove unnecessary data to further reduce the object size
ec_rosmap$RNA$counts <- NULL
ec_rosmap$RNA$data <- NULL
ec_rosmap$RNA$scale.data <- Matrix::Matrix(
    0,
    nrow = nrow(ec_rosmap$RNA),
    ncol = ncol(ec_rosmap$RNA)
)

ec_rosmap@meta.data <-
    ec_rosmap@meta.data |>
    tibble::rownames_to_column("barcode") |>
    dplyr::select(
        barcode,
        rosmap_label,
        rosmap_score
    ) |>
    tibble::column_to_rownames("barcode")

ec_rosmap@commands <- list()
ec_rosmap@tools <- list()

# store colors in Seurat object
set.seed(123)
my_cols_50 <- unname(Polychrome::createPalette(50, pals::cols25()))
ec_rosmap@misc$sample_cols <- my_cols_50
str(ec@meta.data)
scMisc::lss()

qsave(ec_rosmap, file.path("docs", "ec_rosmap.qs"))

pns_sn_sciatic_milbrandt_subset <- subset(
    pns_sn_sciatic_milbrandt,
    subset = cluster %in% c("mySC", "nmSC", "PnC")
)
Idents(pns_sn_sciatic_milbrandt_subset) <- factor(
    pns_sn_sciatic_milbrandt_subset$cluster,
    levels = c("mySC", "nmSC", "PnC")
)

dotPlot(
    path = file.path("lookup", "markers.csv"),
    object = pns_sn_sciatic_milbrandt_subset,
    par = "novel",
    dot_min = 0.07,
    height = 2,
    width = 4,
    ortho = "human2mouse"
)

## Supplementary Figure 2D ----
pns_sn_sciatic_milbrandt <- qs::qread(
    "/home/mischko/Documents/beruf/forschung/scRNA_reference/pns_atlas_milbrandt/pns_sn_sciatic_GSE182098.qs"
)
pns_sn_sciatic_milbrandt_subset <- subset(
    pns_sn_sciatic_milbrandt,
    subset = cluster %in% c("mySC", "nmSC", "PnC")
)

Idents(pns_sn_sciatic_milbrandt_subset) <- factor(
    pns_sn_sciatic_milbrandt_subset$cluster,
    levels = c("mySC", "nmSC", "PnC")
)

scMisc::dotPlot(
    path = file.path("lookup", "markers.csv"),
    object = pns_sn_sciatic_milbrandt_subset,
    par = "novel",
    dot_min = 0.07,
    height = 2,
    width = 4,
    ortho = "human2mouse"
)

dotplot_data_milbrandt <-
    DotPlotData(
        object = pns_sn_sciatic_milbrandt_subset,
        features = c("Mlip", "Grik3", "Prima1", "Cxcl14"),
        dot.min = 0.07,
    )

dotplot_data_milbrandt <-
    DotPlotData(
        object = sc_merge,
        features = c("Mlip", "Grik3", "Prima1", "Cxcl14"),
        dot.min = 0.07,
    )

DotPlotModified(
    data.plot = dotplot_data_milbrandt,
    scale.by = "size",
    dot.scale = 10,
) +
    viridis::scale_color_viridis(option = "viridis") +
    theme(
        axis.text.x = element_text(
            angle = 90,
            vjust = 0.5,
            hjust = 1,
            face = "italic",
            size = 7
        ),
        legend.position = "top",
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 10)
    ) +
    xlab("") +
    ylab("")

# subset sc_merge
sc_merge_subset <- subset(
    sc_merge,
    subset = cluster %in% c("mySC", "nmSC", "periC1", "periC2", "periC3")
)

# rename periC1, periC2, periC3 to periC
sc_merge_subset$cluster <- gsub(
    pattern = "periC\\d",
    replacement = "periC",
    x = sc_merge_subset$cluster
)
sc_merge_subset$cluster <- factor(
    sc_merge_subset$cluster,
    levels = c("mySC", "nmSC", "periC")
)
Idents(sc_merge_subset) <- sc_merge_subset$cluster

dotplot_data_heming <-
    DotPlotData(
        object = sc_merge_subset,
        features = c("MLIP", "GRIK3", "PRIMA1", "CXCL14"),
        dot.min = 0.07,
    )

dotplot_human_rodent <-
    list(
        rodent = dotplot_data_milbrandt,
        human = dotplot_data_heming
    )

qsave(dotplot_human_rodent, file.path("docs", "dotplot_human_rodent.qs"))

# Supplementary Figure 4 ----
## Supplementary Figure 4A ----
dotplot_data_ic <-
    DotPlotData(
        object = ic,
        features = ic_figure@misc$markers_dotplot_ic,
        dot.min = 0.01,
    )

qsave(dotplot_data_ic, file.path("docs", "dotplot_data_ic.qs"))

# Supplementary Figure 4B ----
enrichr_de_clusters <- c("mySC", "nmSC", "PC2")
names(enrichr_de_clusters) <- c("mySC", "nmSC", "PC2")

enrichr_de_pos <-
    lapply(
        enrichr_de_clusters,
        function(condition) {
            readxl::read_excel(
                file.path(
                    "results",
                    "enrichr",
                    paste0("enrichr_pos_", condition, ".xlsx")
                ),
                sheet = "GO_Biological_Process_2023"
            )
        }
    )

enrichr_de_neg <-
    lapply(
        enrichr_de_clusters,
        function(condition) {
            readxl::read_excel(
                file.path(
                    "results",
                    "enrichr",
                    paste0("enrichr_neg_", condition, ".xlsx")
                ),
                sheet = "GO_Biological_Process_2023"
            )
        }
    )
glimpse(enrichr_de)
enrichr_de <-
    list(
        pos = enrichr_de_pos,
        neg = enrichr_de_neg
    )

qsave(enrichr_de, file.path("docs", "enrichr_de.qs"))

# Supplementary Figure 5 ----
# Load and join g-ratio measurements with sample metadata ----
g_ratio <-
    read_excel(file.path("lookup", "g_ratio.xlsx")) |>
    left_join(sample_lookup, join_by(sample)) |>
    mutate(level2 = factor(level2, levels = umap_figure@misc$level2_order)) |>
    select(sample, level2, g_ratio, axon_diameter, ncv_tibial_motoric)

qsave(g_ratio, file.path("docs", "g_ratio.qs"))

## Supplementary Figure 5D ----
axon_count_table <-
    readxl::read_xlsx(file.path("lookup", "axon_count_v2.xlsx")) |>
    dplyr::filter(is.na(remove)) |>
    dplyr::filter(sample != "NA") |>
    mutate(across(
        c(fascicle, normal_myelin:total_myelinated_axons, area_micrometer),
        parse_number
    )) |>
    mutate(fascicle = str_extract(fascicle, "[0-9]+"))

axon_count_fascicle <-
    axon_count_table |>
    group_by(sample, fascicle) |>
    summarize(
        across(c(normal_myelin:total_myelinated_axons, area_micrometer), sum),
        .groups = "drop"
    ) |>
    left_join(sample_lookup, by = "sample") |>
    mutate(level2 = factor(level2, levels = umap_figure@misc$level2_order)) |>
    group_by(sample) |>
    mutate(
        axon_count = mean(normal_myelin),
        ncv_tibial_motoric = mean(ncv_tibial_motoric)
    ) |>
    select(sample, level2, axon_count, ncv_tibial_motoric) |>
    mutate(log_axon_normal = log(axon_count)) |>
    distinct()

qsave(axon_count_fascicle, file.path("docs", "axon_count.qs"))

## Supplementary Figure 5G ----
abundance_axon <-
    table(sc_merge$cluster, sc_merge$sample) |>
    as.data.frame.matrix() |>
    rownames_to_column("cell") |>
    mutate(across(where(is.numeric), function(x) x / sum(x) * 100)) |>
    pivot_longer(!cell, names_to = "sample", values_to = "count") |>
    left_join(axon_count_fascicle, join_by(sample))

qsave(abundance_axon, file.path("docs", "abundance_axon.qs"))

## Supplementary Figure 5I ----
xenium_objects <- qs::qread(file.path("objects", "xenium_objects.qs"))

cells_list <-
    lapply(xenium_objects, function(x) (x$sn_predictions)) |>
    unlist()

cells_predicted <-
    tibble(cluster = unname(cells_list), sample = names(cells_list)) |>
    mutate(condition = str_replace(sample, "(S\\d+)_(\\w+).*", "\\2")) |>
    mutate(sample = str_replace(sample, "(S\\d+)_.*", "\\1")) |>
    dplyr::filter(!is.na(cluster))


# boxplot sc
xenium_sc_t_nk <-
    cells_predicted |>
    dplyr::count(cluster, sample) |>
    pivot_wider(names_from = sample, values_from = n) |>
    mutate(across(
        where(is.numeric),
        function(x) x / sum(x, na.rm = TRUE) * 100
    )) |>
    pivot_longer(!cluster, names_to = "sample", values_to = "percent") |>
    left_join(select(sample_lookup, sample, level2)) |>
    dplyr::rename(condition = level2) |>
    dplyr::filter(cluster %in% c("mySC", "nmSC", "repairSC", "T_NK")) |>
    mutate(percent = replace_na(percent, 0)) |>
    mutate(condition = factor(condition, levels = sc_merge@misc$level2_order))

qsave(xenium_sc_t_nk, file.path("docs", "xenium_sc_t_nk.qs"))
