#===============================================================================
# Differential Abundance Analysis using PCA
#===============================================================================
# Purpose: Perform Principal Component Analysis (PCA) on cluster abundances in 
# single-cell data to identify patterns and visualize differences between samples.
#
# Methods: 
# - Load preprocessed single-cell data
# - Define a function to compute and plot PCA results
# - Apply the function to specific subsets of data
#===============================================================================

# Libraries ----
# Load necessary libraries for data manipulation and visualization
library(Seurat)
library(tidyverse)
library(factoextra)
library(FactoMineR)
library(viridis)
library(ggrepel)

# General settings ----
# Resolve function conflicts
conflicts_prefer(base::setdiff)
conflicts_prefer(base::as.data.frame)

# Load preprocessed data ----
# Read in the Seurat objects containing single-cell data
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"))
ic <- qs::qread(file.path("objects", "ic.qs"))

table(ic$level2)
table(ic_cidp_vn_ctrl_ciap$level2)

# Define PCA plotting function ----
# Function to perform PCA on cluster abundances and generate plots
pcaSeurat <- function(object, label1, label2, label3, label3_colors = NULL) {
  object_parse <- deparse(substitute(object))
  
  # Create a matrix of cluster abundances per sample
  cl_size <-
    base::as.data.frame.matrix(table(object@meta.data[[label1]], object@meta.data[[label2]])) |>
    t()
  
  colnames(cl_size) <- levels(object@meta.data[[label1]])
  
  # Perform PCA on the abundance matrix
  pca_result <- FactoMineR::PCA(cl_size, scale.unit = TRUE, ncp = 30, graph = FALSE)
  
  # Plot the variance explained by each principal component
  factoextra::fviz_eig(pca_result, addlabels = TRUE, ylim = c(0,50), ncp = 7)
  ggsave(paste0("./results/pca/pca_", object_parse, "_", label3, "_eigen.pdf"))
  
  # Plot the variables contributing to the principal components
  pca_var_plot <-
    factoextra::fviz_pca_var(
      pca_result,
      col.var = "contrib",
      gradient.cols = viridis::plasma(100),
      repel = TRUE,
      select.var = list(contrib = 10)
    ) +
    labs(
      title = "",
      x = "PC1", y = "PC2"
    ) +
    theme_classic() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_rect(color = "black", size = 1, fill = NA),
      aspect.ratio = 1
    )
  
  # Prepare metadata for plotting
  lookup_pre <-
    data.frame(
      label1 = object@meta.data[[label2]],
      label3 = object@meta.data[[label3]]
    ) |>
    distinct()
  
  lookup <-
    data.frame(label1 = rownames(cl_size)) |>
    left_join(lookup_pre)
  
  # Plot PCA results with samples colored by a specified label
  pca_plot_group <-
    factoextra::fviz_pca_ind(
      pca_result,
      pointsize = 5,
      pointshape = 21,
      geom.ind = "point",
      fill.ind = lookup$label3,
      col.ind = "black",
      palette = label3_colors,
      addEllipses = FALSE,
      ellipse.type = "confidence",
      legend.title = label3,
      axes.linetype = "solid",
      alpha.ind = 1
    )
  
  pca_ggplot_group <-
    ggpubr::ggpar(
      pca_plot_group,
      title = "",
      xlab = "PC1",
      ylab = "PC2",
      ticks = FALSE,
      tickslab = FALSE,
      ggtheme = theme_classic() +
        theme(
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          plot.title = element_text(size = 25),
          panel.border = element_rect(color = "black", size = 1, fill = NA),
          aspect.ratio = 1,
        )
    ) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
    ) +
    geom_text(aes(label = lookup$label1), nudge_y = .3)
  
  if(is.numeric(lookup$label3)) {
    pca_ggplot_group <- pca_ggplot_group + viridis::scale_fill_viridis(option = "magma")
  }
  
  # Combine variable and individual plots
  pca_plots <- patchwork::wrap_plots(pca_var_plot, pca_ggplot_group, ncol = 2)
  
  # Save the combined plots to a file
  ggsave(paste0("./results/pca/pca_", object_parse, "_", label3,  ".pdf"), width = 12, height = 6,
         plot = pca_plots)
}

debugonce(pcaSeurat)

# Run PCA abundance function ----
# Subset data to include specific conditions
vn_cidp_ciap_ctrl <- subset(sc_merge, level2 %in% c("VN", "CIDP", "CIAP", "CTRL"))

# Perform PCA and plot results for different labels
pcaSeurat(
  object = vn_cidp_ciap_ctrl,
  label1 = "cluster",
  label2 = "sample",
  label3 = "level2",
  label3_colors = vn_cidp_ciap_ctrl@misc$level2_cols
)

pcaSeurat(
  object = vn_cidp_ciap_ctrl,
  label1 = "cluster",
  label2 = "sample",
  label3 = "log_axon_normal",
)

pcaSeurat(
  object = vn_cidp_ciap_ctrl,
  label1 = "cluster",
  label2 = "sample",
  label3 = "g_ratio"
)

pcaSeurat(
  object = vn_cidp_ciap_ctrl,
  label1 = "cluster",
  label2 = "sample",
  label3 = "axon_diameter"
)

# Convert 'incat' to numeric for analysis
vn_cidp_ciap_ctrl$incat <- as.numeric(vn_cidp_ciap_ctrl$incat)

pcaSeurat(
  object = vn_cidp_ciap_ctrl,
  label1 = "cluster",
  label2 = "sample",
  label3 = "incat"
)

pcaSeurat(
  object = vn_cidp_ciap_ctrl,
  label1 = "cluster",
  label2 = "sample",
  label3 = "mrvi_cluster",
  label3_colors =  pals::cols25()
)
