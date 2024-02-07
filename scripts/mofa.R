# libraries  ----
library(Seurat)
library(BPCells)
library(SeuratObject)
library(tidyverse)
library(qs)
library(MOFAcellulaR)
library(conflicted)

# general settings  ----
options(warn = 0)
options(Seurat.object.assay.version = "v5")
future::plan("multicore", workers = 6)
my_cols_25 <- pals::cols25()
my_cols_50 <- unname(Polychrome::createPalette(50, pals::cols25()))
conflicts_prefer(base::setdiff)
conflicts_prefer(base::unname)

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)

sc_merge$cluster_sample <- paste0(sc_merge$cluster, "_", sc_merge$sample)

# pseudobulk expression 
bulk <- AverageExpression(
    sc_merge,
    method = "aggregate",
    return.seurat = FALSE,
    slot = "counts",
    assays = "RNA",
    group.by = "cluster_sample"
)

# _ was replaced by - 
colnames(bulk$RNA) <- as.character(colnames(bulk$RNA))
colnames(bulk$RNA) <- gsub("-", "_", colnames(bulk$RNA))


cluster_sample_bulk <- as.character(colnames(bulk$RNA))
count_cells <- dplyr::count(sc_merge@meta.data, cluster_sample)

# get metadata from Seurat object and preprocess
# be careful with the order of the metadata
metadata_mofa <- sc_merge@meta.data |>
    select(cluster_sample, sample, cluster, level1, level2, sex, age, center, incat) |>
    distinct(cluster_sample, .keep_all = TRUE) |>
    left_join(count_cells, by = "cluster_sample") |>
    arrange(match(cluster_sample, cluster_sample_bulk)) |>
    as_tibble()

cluster_of_interest <- unique(as.character(metadata_mofa$cluster))


# sanity check
all.equal(metadata_mofa$cluster_sample, cluster_sample_bulk)

# # mofa analysis
# following this tutorial https://saezlab.github.io/MOFAcellulaR/articles/get-started.html#exporting-model-outputs
inputs_dir <- base::system.file("extdata", package = "MOFAcellulaR")
load(file.path(inputs_dir, "testpbcounts.rda"))
load(file.path(inputs_dir, "testcoldata.rda"))

testcoldata %>%
    dplyr::select(donor_id) %>%
    dplyr::group_by(donor_id) %>%
    dplyr::summarise(n()) %>%
    head()

# testpbcounts[1:5, 1:5]
# processing pseudobulk expression ----
pb_obj <- MOFAcellulaR::create_init_exp(counts = testpbcounts, coldata = testcoldata)

pb_obj <- MOFAcellulaR::create_init_exp(counts = bulk$RNA, coldata = metadata_mofa)

ct_list <- MOFAcellulaR::filt_profiles(
    pb_dat = pb_obj,
    cts = c("mySC", "nmSC", "repairSC", "damageSC"),
    # cts = cluster_of_interest,
    ncells = 10, # Change to your knowledge!!
    counts_col = "n", # This refers to the column name in testcoldata where the number of cells per profile was stored
    ct_col = "cluster"
) # This refers to the column name in testcoldata where the cell-type label was stored


ct_list_backup <- MOFAcellulaR::filt_profiles(pb_dat = pb_obj,
                          cts = c("Fib","CM"),
                          ncells = 0, # Change to your knowledge!! 
                          counts_col = "cell_counts", # This refers to the column name in testcoldata where the number of cells per profile was stored
                          ct_col = "cell_type") # This refers to the column name in testcoldata where the cell-type label was stored

#sanity check
colnames(ct_list$Macro1)

# filter views with few samples
ct_list <- MOFAcellulaR::filt_views_bysamples(
    pb_dat_list = ct_list,
    nsamples = 2
)

# identify lowly expressed genes
ct_list <- MOFAcellulaR::filt_gex_byexpr(
    pb_dat_list = ct_list,
    min.count = 5, # Modify!!
    min.prop = 0.25
) # Modify!!


# filter view with few genes
ct_list <- filt_views_bygenes(
    pb_dat_list = ct_list,
    ngenes = 10
)


# somehow this only works whenn function is redefined

my_filt_samples_bycov <- function(pb_dat_list, prop_coverage = 0.9) {
    pb_dat_list_new <- purrr::map(pb_dat_list, function(pb_obj) {
        mat <- SummarizedExperiment::assay(pb_obj, "counts")
        sample_vect <- (mat != 0) %>% colSums(.) / nrow(mat)
        sel_samples <- names(sample_vect[which(sample_vect >= prop_coverage)])
        return(pb_obj[, sel_samples])
    })
    return(pb_dat_list_new)
}

# filter samples with low coverage
ct_list <- my_filt_samples_bycov(
    pb_dat_list = ct_list,
    prop_coverage = 0.9
)

# normalization via TMM
ct_list <- MOFAcellulaR::tmm_trns(
    pb_dat_list = ct_list,
    scale_factor = 1000000
)

# identify highly variable genes
ct_list <- MOFAcellulaR::filt_gex_byhvg(
    pb_dat_list = ct_list,
    prior_hvg = NULL,
    var.threshold = 0
)

# prior_hvg_test <- list(
#     "CM" = c("TTN"),
#     "Fib" = c("POSTN")
# )

# ct_list <- MOFAcellulaR::filt_gex_bybckgrnd(
#     pb_dat_list = ct_list,
#     prior_mrks = prior_hvg_test
# )

# filter views with few genes
ct_list <- MOFAcellulaR::filt_views_bygenes(
    pb_dat_list = ct_list,
    ngenes = 15
)

# convert into MOFA object
multiview_dat <- pb_dat2MOFA(
    pb_dat_list = ct_list,
    sample_column = "sample"
)

# fitting a MOFA model ---
MOFAobject <- MOFA2::create_mofa(multiview_dat)

data_opts <- MOFA2::get_default_data_options(MOFAobject)
train_opts <- MOFA2::get_default_training_options(MOFAobject)
model_opts <- MOFA2::get_default_model_options(MOFAobject)

# This avoids the regularization of multicellular programs per cell type.
# This avoids less sparse gene weights
model_opts$spikeslab_weights <- FALSE

# Define the number of factors needed
model_opts$num_factors <- 5

# Prepare MOFA model:
MOFAobject <- MOFA2::prepare_mofa(
    object = MOFAobject,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
)

outfile <- file.path("objects", "MOFAobject.rds")
outfile <- file.path("objects", "MOFAobject_SC.rds")

# reticulate::use_condaenv()
reticulate::py_config()

model <- MOFA2::run_mofa(MOFAobject, outfile, use_basilisk = TRUE)

# # exploring the MOFA model ----
# metadata <- readRDS(file.path(inputs_dir, "testmetadata.rds"))

# set.seed(145)
# metadata$fake_var <- stats::rnorm(nrow(metadata))

all_factors <- MOFAcellulaR::get_tidy_factors(
    model = model,
    metadata = metadata_mofa,
    factor = "all",
    sample_id_column = "sample"
)

# umap
umap_embedding <- MOFAcellulaR::plot_sample_2D(
    model = model,
    method = "UMAP",
    metadata = metadata_mofa_sample,
    sample_id_column = "sample",
    color_by = "level2"
)

umap_plot <-
    ggplot2::ggplot(umap_embedding, aes(x = UMAP_1, y = UMAP_2, color = color_col)) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 12)) +
    ggplot2::labs(color = "") + 
    scale_color_manual(values = my_cols_50) + 
    geom_text(aes(label = sample), vjust = 2) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        aspect.ratio = 1,
    )

ggsave(plot = umap_plot, file = file.path("results", "mofa", "mofa_umap.pdf"), width = 10, height = 8)
ggsave(plot = umap_plot, file = file.path("results", "mofa", "mofa_sc_umap.pdf"), width = 10, height = 8)

metadata_mofa_sample <- sc_merge@meta.data |>
    select(cluster_sample, sample, cluster, level1, level2, sex, age, center, incat) |>
    distinct(sample, .keep_all = TRUE) |>
    left_join(count_cells, by = "cluster_sample") |>
    arrange(match(cluster_sample, cluster_sample_bulk)) |>
    mutate(level0 = if_else(level2 == "CTRL", "CTRL", "PNP")) |>
    as_tibble()

dplyr::count(metadata_mofa_sample, level0, level2)

# mds
mds_embedding <- MOFAcellulaR::plot_sample_2D(
    model = model,
    method = "MDS",
    metadata = metadata_mofa_sample,
    sample_id_column = "sample",
    color_by = "level2"
)

mds_plot <-
    ggplot2::ggplot(mds_embedding, aes(x = MDS1, y = MDS2, color = color_col)) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 12)) +
    ggplot2::labs(color = "") + 
    scale_color_manual(values = my_cols_50) + 
    geom_text(aes(label = sample), vjust = 2) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        aspect.ratio = 1,
    )

ggsave(plot = mds_plot, file = file.path("results", "mofa", "mofa_mds.pdf"), width = 10, height = 8)

categorical_assoc <- MOFAcellulaR::get_associations(
    model = model,
    metadata = metadata_mofa_sample,
    sample_id_column = "sample",
    test_variable = "level0",
    test_type = "categorical",
    group = FALSE
)

continuous_assoc <- MOFAcellulaR::get_associations(
    model = model,
    metadata = metadata_mofa_sample,
    sample_id_column = "sample",
    test_variable = "incat",
    test_type = "categorical",
    group = FALSE
)

assoc_list <- list("disease" = categorical_assoc, "incat" = continuous_assoc)

plot_MOFA_hmap(
    model = model,
    group = FALSE,
    metadata = metadata_mofa_sample,
    sample_id_column = "sample",
    sample_anns = c("level0", "level2", "sample"),
    assoc_list = assoc_list
)

ggsave(file = file.path("results", "mofa", "mofa_hmap.pdf"), width = 10, height = 8)
