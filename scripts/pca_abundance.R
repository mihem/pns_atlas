#libraries ---
library(Seurat)
library(tidyverse)
library(factoextra)
library(FactoMineR)
library(viridis)
library(ggrepel)

# general settings ----
conflicts_prefer(base::setdiff)
conflicts_prefer(base::as.data.frame)

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"))
ic <- qs::qread(file.path("objects", "ic.qs"))

table(ic$level2)
table(ic_cidp_vn_ctrl_ciap$level2)

# function to plot pca of Seurat abundancies
pcaSeurat <- function(object, label1, label2, label3, label3_colors = NULL) {
  object_parse <- deparse(substitute(object))
  cl_size <-
    base::as.data.frame.matrix(table(object@meta.data[[label1]], object@meta.data[[label2]])) |>
    t()

  colnames(cl_size) <- levels(object@meta.data[[label1]])

  pca_result <- FactoMineR::PCA(cl_size, scale.unit = TRUE, ncp = 30, graph = FALSE)
  factoextra::fviz_eig(pca_result, addlabels = TRUE, ylim = c(0,50), ncp = 7)
  ggsave(paste0("./results/pca/pca_", object_parse, "_", label3, "_eigen.pdf"))

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

  lookup_pre <-
    data.frame(
      label1 = object@meta.data[[label2]],
      label3 = object@meta.data[[label3]]
    ) |>
    distinct()

  lookup <-
    data.frame(label1 = rownames(cl_size)) |>
    left_join(lookup_pre)

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

  pca_plots <- patchwork::wrap_plots(pca_var_plot, pca_ggplot_group, ncol = 2)
  ggsave(paste0("./results/pca/pca_", object_parse, "_", label3,  ".pdf"), width = 12, height = 6,
         plot = pca_plots)

}

debugonce(pcaSeurat)

# run pca abundance function ----
vn_cidp_ciap_ctrl <- subset(sc_merge, level2 %in% c("VN", "CIDP", "CIAP", "CTRL"))

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

# ic_vn_cidp_ciap_ctrl <- subset(ic, level2 %in% c("CIDP", "VN", "CTRL", "CIAP"))

# pcaSeurat(
#   object = ic_vn_cidp_ciap_ctrl,
#   label1 = "ic_cluster",
#   label2 = "sample",
#   label3 = "level2",
#   label3_colors = ic@misc$level2_cols
# )

# pcaSeurat(
#   object = ic_vn_cidp_ciap_ctrl,
#   label1 = "ic_cluster",
#   label2 = "sample",
#   label3 = "log_axon_normal",
# )

# pcaSeurat(
#   object = ic_vn_cidp_ciap_ctrl,
#   label1 = "ic_cluster",
#   label2 = "sample",
#   label3 = "g_ratio"
# )

# pcaSeurat(
#   object = ic_vn_cidp_ciap_ctrl,
#   label1 = "ic_cluster",
#   label2 = "sample",
#   label3 = "axon_diameter"
# )

# ic_vn_cidp_ciap_ctrl$incat <- as.numeric(ic_vn_cidp_ciap_ctrl$incat)

# pcaSeurat(
#   object = ic_vn_cidp_ciap_ctrl,
#   label1 = "ic_cluster",
#   label2 = "sample",
#   label3 = "incat"
# )

# NMF not superior to PCA
# library(NMF)

# cl_size_sc <-
#   as.data.frame.matrix(table(sc_merge@meta.data[["cluster"]], sc_merge@meta.data[["sample"]])) |>
#   t()
# colnames(cl_size_sc) <- levels(sc_merge@meta.data[["cluster"]])

# nmf_sc <- NMF::nmf(cl_size_sc, 2, method = "Frobenius")

# coldata <- readr::read_csv(file.path("lookup", "sample_lookup.csv")) 

# nmf_df <-
#  as.data.frame(NMF::basis(nmf_sc)) |>
#  mutate(level2 = coldata$level2) |>
#  tibble::rownames_to_column("sample")


# ggplot(nmf_df, aes(x = V1, V2, color = sample)) +
#   geom_point()
