# plot quality control metrics

# libraries ---
library(tidyverse)
library(Seurat)
library(BPCells)
library(pals)

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)

my_cols_50 <- unname(Polychrome::createPalette(50, pals::cols25()))

count_cells <- dplyr::count(sc_merge@meta.data, sample)

count_genes <-
    dplyr::tibble(feature = sc_merge@meta.data$nFeature_RNA, sample = sc_merge@meta.data$sample) |>
    dplyr::group_by(sample) |>
    dplyr::summarize(median_genes = median(feature))


count_genes <- dplyr::tibble(
    feature = sc_merge@meta.data$nFeature_RNA,
    sample = sc_merge@meta.data$sample
)

count_genes_plot <-
    count_genes |>
    ggplot(aes(x = sample, y = feature, fill = sample)) +
    geom_boxplot() + 
    scale_fill_manual(values = my_cols_50) + 
    theme_bw() + 
    theme(legend.position = "none") + 
    xlab("") + 
    ylab("") + 
    ggtitle("Number of genes per nucleus")

ggsave(
    file.path("results", "qc", "count_genes.pdf"),
    plot = count_genes_plot,
    width = 10,
    height = 5
)
