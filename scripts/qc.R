#===============================================================================
# Quality Control Metrics Plotting
#===============================================================================
# Purpose: Generate and save plots for quality control metrics of single-cell RNA 
# sequencing data.
#===============================================================================

# Load necessary libraries
library(tidyverse)
library(Seurat)
library(BPCells)
library(pals)

# Load preprocessed data
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)

# Define color palette
my_cols_50 <- unname(Polychrome::createPalette(50, pals::cols25()))

# Count cells per sample
count_cells <- dplyr::count(sc_merge@meta.data, sample)

# Calculate median number of genes per sample
count_genes <- dplyr::tibble(
    feature = sc_merge@meta.data$nFeature_RNA,
    sample = sc_merge@meta.data$sample
) |>
dplyr::group_by(sample) |>
dplyr::summarize(median_genes = median(feature))

# Create a boxplot for the number of genes per nucleus
count_genes_plot <- count_genes |>
    ggplot(aes(x = sample, y = feature, fill = sample)) +
    geom_boxplot() + 
    scale_fill_manual(values = my_cols_50) + 
    theme_bw() + 
    theme(legend.position = "none") + 
    xlab("") + 
    ylab("") + 
    ggtitle("Number of genes per nucleus")

# Save the plot
ggsave(
    file.path("results", "qc", "count_genes.pdf"),
    plot = count_genes_plot,
    width = 10,
    height = 5
)
