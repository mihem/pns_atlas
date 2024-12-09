#===============================================================================
# Preprocessing Single-Cell Data
#===============================================================================
# Purpose: Preprocess single-cell RNA sequencing data for downstream analysis.
#
# Methods: 
# - Load raw data
# - Perform quality control and filtering
# - Normalize and scale data
# - Identify highly variable features
# - Perform dimensionality reduction
# - Save the preprocessed data for further analysis
#===============================================================================

# libraries required for analysis ----
library(Seurat)
library(SeuratWrappers)
library(BPCells)
library(Azimuth)
library(SeuratObject)
library(tidyverse)
library(writexl)
library(patchwork)
library(conflicted)
library(scDblFinder)
library(qs)
library(scMisc)
library(janitor)
library(biomaRt)
library(readxl)

# optional libraries for better coding
library(httpgd)
library(languageserver)
library(codegrip)

install.packages("languageserver")

# general settings  ----
options(warn = 0)
options(Seurat.object.assay.version = "v5")
future::plan("multicore", workers = 6)
options(future.globals.maxSize = 16000 * 1024^2)
conflicts_prefer(base::setdiff)

# create folders ---
folders <- c(
  "qc",
  "table",
  "umap",
  "dotplot",
  "abundance",
  "de",
  "enrichr",
  "map",
  "pca",
  "miloDE",
  "miloR",
  "cna",
  "mofa",
  "histo",
  "module",
  "projectil",
  "venn",
  "gratio",
  "demographics",
  "violinplot"
)
lapply(file.path("results", folders), dir.create, recursive = TRUE)

# meta data ----
#calculcate age or use lookup
sample_lookup <- 
  readr::read_csv(file.path("lookup", "sample_lookup.csv")) |>
  janitor::clean_names() |>
  mutate(age_calc = lubridate::time_length(difftime(nerve_date, birth_date), "years")) |>
  mutate(age_calc = floor(age_calc)) |>
  mutate(age = coalesce(age_calc, age)) |>
  dplyr::select(-age_calc)

axon_counts <- 
  readxl::read_excel(file.path("lookup", "axon_count_v2.xlsx"), na = c("", "NA"))

# create overview table
overview_table <-
    sample_lookup |>
    dplyr::mutate(sex_cat = if_else(sex == "male", 1, 0)) |>
    dplyr::group_by(level2) |>
    dplyr::summarize(n = n(),
                     age = mean(age, na.rm = TRUE),
                     female = (1 - mean(sex_cat, na.rm = TRUE))*100
                     )
writexl::write_xlsx(overview_table, file.path("results", "table", "overview_table.xlsx"))

incat <- 
  sample_lookup |>
  select(sample, level2, incat) |>
  dplyr::group_by(level2) |>
  dplyr::filter(level2 != "CTRL") |>
  dplyr::mutate(incat = as.numeric(incat)) |>
  summarize(mean = mean(incat, na.rm = TRUE))

writexl::write_xlsx(incat, file.path("results", "table", "incat.xlsx"))

# find files and match to names from lookup table
h5_path <- list.files(pattern = ".h5", recursive = TRUE)

sample_names <-
  tibble(internal_name = h5_path) |>
  dplyr::mutate(internal_name = gsub(x = internal_name, pattern = "([^/]+)/([^/]+)/(.+)", replacement = "\\2")) |>
  left_join(sample_lookup)
    
#sanity check
setdiff(sample_names$internal_name, sample_lookup$internal_name)
setdiff(sample_lookup$internal_name, sample_names$internal_name)

# read in data and create Seurat object -----
sc_raw_list <- lapply(h5_path, BPCells::open_matrix_10x_hdf5) |>
    setNames(sample_names$sample)

# write the matrix to a directory
map2(
  .x = sc_raw_list,
  .y = file.path("matrix", names(sc_raw_list)),
  .f = BPCells::write_matrix_dir
)

# load matrix from directory
sc_raw_mat <-
  map(
    .x = file.path("matrix", names(sc_raw_list)),
    .f = open_matrix_dir
  ) |>
  setNames(sample_names$sample)

# hacky solution because biomaRt timeout
MyConvertEnsembleToSymbol <- function(mat, gene_IDs, symbol = "hgnc_symbol") {
  name_df <- data.frame(gene_id = c(rownames(mat)))
  name_df$orig.id <- name_df$gene_id
  name_df$gene_id <- as.character(name_df$gene_id)
  name_df$gene_id <- sub("[.][0-9]*", "", name_df$gene_id)
  gene.df <- dplyr::left_join(name_df, gene_IDs, by = c(gene_id = "ensembl_gene_id"))
  rownames(gene.df) <- make.unique(gene.df$orig.id)
  gene.df <- gene.df[rownames(mat), ]
  gene.df <- gene.df[gene.df[, symbol] != "", ]
  gene.df <- gene.df[!is.na(gene.df$orig.id), ]
  mat.filter <- mat[gene.df$orig.id, ]
  rownames(mat.filter) <- make.unique(gene.df[, symbol])
  return(mat.filter)
}

# get ensemble genes
genes_human <- rownames(sc_raw_mat[[1]])

# create the biomaRt object
mart_human <- biomaRt::useDataset(dataset = "hsapiens_gene_ensembl", biomaRt::useMart("ensembl"))

# Get the gene IDs and HGNC symbols
gene_IDs_human <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  values = genes_human,
  mart = mart_human
)

# Convert the ensemble gene IDs to HGNC symbols
sc_raw_mat <-
  map(
    .x = sc_raw_mat,
    .f = function(x, y) MyConvertEnsembleToSymbol(mat = x, gene_IDs = gene_IDs_human)
  )

#sanity check
str(rownames(sc_raw_mat$S01), max.level = 2)

sc_list <- map(sc_raw_mat, CreateSeuratObject, min.cells = 3, min.features = 200, project = "sural")

# reorder sc_list based on the sample names
sc_list <- sc_list[order(names(sc_list))]

# plot QC: MT and genes ----
for (i in seq_along(sc_list)) {
    sc_list[[i]][["percent_mt"]] <- PercentageFeatureSet(sc_list[[i]], pattern = "^MT")
}

plot1 <- vector("list", length = length(sc_list))

for (i in seq_along(sc_list)) {
  plot1[[i]] <- FeatureScatter(object = sc_list[[i]], feature1 = "nCount_RNA", feature2 = "percent_mt", raster =  FALSE) +
    labs(title = names(sc_list)[[i]]) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
    NoLegend()
}

plot1_patch <- patchwork::wrap_plots(plot1, ncol = 4)
ggsave(file.path("results", "qc", "mt.png"), width = 15, height = 40, limitsize = FALSE)

plot2 <- vector("list", length = length(sc_list))

for (i in seq_along(sc_list)) {
  plot2[[i]] <- FeatureScatter(object = sc_list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    labs(title = names(sc_list)[[i]]) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
    NoLegend()
}

plot2_patch <- patchwork::wrap_plots(plot2, ncol = 4)
ggsave(file.path("results", "qc", "genes.png"), width = 15, height = 40, limitsize = FALSE)

qs::qsave(sc_list, file.path("objects", "sc_list.qs"))

# doublet detecting with scDblFinder ----

#convert in-memory-matrix with as
#calculcate doubles and add scDblFinder scores to seurat objects
doubletFun <- function(seu_obj) {
  sc_list_on_disk <- as(object = seu_obj[["RNA"]]$counts, Class = "dgCMatrix")
  sce <- scDblFinder(sc_list_on_disk)
  seu_obj$scDblFinder.score <- sce$scDblFinder.score
  seu_obj$scDblFinder.class <- sce$scDblFinder.class
  return(seu_obj)
}

sc_list <- purrr::map(sc_list, doubletFun)

#plot scDblFinder
plot3 <- vector("list")

for (i in seq_along(sc_list)) {
  plot3[[i]] <- FeatureScatter(object = sc_list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "scDblFinder.class") +
    labs(title = names(sc_list)[i]) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))
}

plot3_patch <- patchwork::wrap_plots(plot3, ncol = 4)
ggsave(file.path("results", "qc", "doublet.png"), width = 15, height = 40, limitsize = FALSE)

#create table with doublet rate
doublet_tbl <- purrr::map_dfr(purrr::map(sc_list, "meta.data"), dplyr::count, scDblFinder.class) |>
    dplyr::mutate(sample = rep(names(sc_list), each = 2)) |>
    tidyr::pivot_wider(names_from = "scDblFinder.class", values_from = "n") |>
    dplyr::mutate(doublet_pct = doublet/(singlet+doublet)*100)

write_csv(doublet_tbl, file.path("results", "qc", "doublets.csv"))

# filter low quality cells and doublets ---
filter_df <- readr::read_csv(file.path("lookup", "filter_df.csv"))

sc_filter <- vector("list", length(sc_list))

for (i in seq_along(sc_list)) {
  sc_filter[[i]] <-
    subset(
      sc_list[[i]],
      subset = nFeature_RNA > 200 & nFeature_RNA < filter_df$rna[[i]] & percent_mt < filter_df$mt[[i]] & scDblFinder.class == "singlet"
    )
}

names(sc_filter) <- names(sc_list)

# check filter settings second round ---

plot4 <- vector("list", length = length(sc_filter))

for (i in seq_along(sc_filter)) {
  plot4[[i]] <- FeatureScatter(object = sc_filter[[i]], feature1 = "nCount_RNA", feature2 = "percent_mt", raster =  FALSE) +
    labs(title = names(sc_filter)[[i]]) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
    NoLegend()
}

plot4_patch <- patchwork::wrap_plots(plot4, ncol = 4)
ggsave(file.path("results", "qc", "mt_post.png"), width = 15, height = 40, limitsize = FALSE)

plot5 <- vector("list", length = length(sc_filter))

for (i in seq_along(sc_filter)) {
  plot5[[i]] <- FeatureScatter(object = sc_filter[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    labs(title = names(sc_filter)[[i]]) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
    NoLegend()
}

plot5_patch <- patchwork::wrap_plots(plot5, ncol = 4)
ggsave(file.path("results", "qc", "genes_post.png"), width = 15, height = 40, limitsize = FALSE)

# merge seurat objects  ------------------------------------------
sc_merge_pre <- merge(x = sc_filter[[1]], y = sc_filter[-1], merge.data = TRUE, add.cell.ids = names(sc_list))
sc_merge_pre$sample <- str_extract(colnames(sc_merge_pre), pattern = "[^_]+")

#this needs to be rejoined and then split again to get the correct names (better naming and required for integration)
sc_merge_pre <- JoinLayers(sc_merge_pre)
sc_merge_pre <- split(x = sc_merge_pre, f = sc_merge_pre$sample)

# add metadata
sample_lookup_sel <-
  sample_lookup |>
  dplyr::select(sample, level1, level2, sex, age, center, species, disease_duration_in_months, incat, mrc_sum_score_60, csf_protein)

sc_merge@meta.data <-
    sc_merge@meta.data |>
    tibble::rownames_to_column("barcode") |>
    dplyr::left_join(sample_lookup_sel, by = "sample") |>
    tibble::column_to_rownames(var = "barcode")


# qc metrics  -----
metrics_files <- list.files("raw", pattern = "metrics_summary.csv", recursive = TRUE, full.names = TRUE)

metrics_data <-
  purrr::map_df(metrics_files, read_csv) |>
  mutate(sample = sample_names$sample, .before = `Estimated Number of Cells`) |>
  arrange(sample)

write_csv(metrics_data, file.path("results", "qc", "cellranger_metrics.csv"))

count_cells <-
  purrr::map_df(sc_list, dim) |>
  dplyr::slice(2) |>
  tidyr::pivot_longer(everything(), names_to = "sample") |>
  dplyr::left_join(dplyr::count(sc_merge_pre@meta.data, sample)) |>
  dplyr::rename(before = value, after = n)

write_csv(count_cells, file.path("results", "qc", "count_cells.csv"))

count_genes <-
    dplyr::bind_cols(feature = sc_merge_pre@meta.data$nFeature_RNA, sample = sc_merge_pre@meta.data$sample) |>
    dplyr::group_by(sample) |>
    dplyr::summarize(median_genes_after = median(feature))

write_csv(count_genes, file.path("results", "qc", "count_genes.csv"))

# important!, fix issue with nested matrix -> leads to problems with file size, NormalizeData and FindTransferAnchors
# https://github.com/satijalab/seurat/issues/7373

for (i in Layers(sc_merge_pre)){
    sc_merge_pre[["RNA"]][[i]] <- as(sc_merge_pre[["RNA"]][[i]], "dgCMatrix")
    sc_merge_pre[["RNA"]][[i]] <- write_matrix_dir(sc_merge_pre[["RNA"]][[i]], dir = paste0("matrix_final/", gsub("counts.", "", i)))
}

# sanity check
str(sc_merge_pre$RNA$counts@matrix@matrix@matrix_list[[1]]@matrix@matrix_list[[1]]@matrix@matrix, max.level = 2)
str(sc_merge_pre$RNA$counts@matrix@matrix, max.level = 2)

# save data
qs::qsave(sc_merge_pre, file.path("objects", "sc_merge_pre.qs"))
